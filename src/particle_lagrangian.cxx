/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <iomanip>
#include <cmath>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "constants.h"

#include <hdf5.h>
#include "hdf5_interface.h"

#include "particle_lagrangian.h"
#include "particle_lagrangian_kernels.h"

namespace plk = Particle_lagrangian_kernels;

namespace
{
    template<typename T>
    void smart_resize(
            std::vector<T>& v,
            const int new_size,
            const double margin)
    {
        const int size     = v.size();      // Used elements.
        const int capacity = v.capacity();  // Reserved size.

        if (new_size * margin < capacity)
        {
            // Too large! Decrease capacity.
            const int new_capacity = int(new_size * margin);
            std::vector<T>(v.begin(), v.begin() + new_size).swap(v);
            v.reserve(new_capacity);
        }
        else if (new_size > capacity)
        {
            // Too small! Increase capacity.
            const int new_capacity = int(new_size * margin);
            v.reserve(new_capacity);
        }

        v.resize(new_size);
    }
}


template<typename TF>
Particle_lagrangian<TF>::Particle_lagrangian(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_particle = inputin.get_item<bool>("particle_lagrangian", "sw_particle", "", false);

    if (sw_particle)
    {
        n_particles = inputin.get_item<int>("particle_lagrangian", "n_particles", "");

        // Raw dump of all particles.
        sw_dump = inputin.get_item<bool>("particle_lagrangian", "sw_dump", "", false);
        if (sw_dump)
        {
            const int sampletime = inputin.get_item<int>("particle_lagrangian", "sampletime_dump", "");
            isampletime_dump = convert_to_itime(sampletime);
        }
    }
}


template<typename TF>
Particle_lagrangian<TF>::~Particle_lagrangian()
{
}


#ifndef USECUDA
template<typename TF>
void Particle_lagrangian<TF>::exec()
{
    // Calculate particle tendencies by tri-linear interpolation of Eulerian velocity fields to particle locations.
    if (!sw_particle)
        return;

    auto& gd = grid.get_grid_data();

    auto diagnose_tendency = [&](
        std::vector<TF>& velocity,
        std::vector<TF>& tendency,
        const std::vector<TF>& fld_3d,
        const std::array<int,3>& loc)
    {
        const TF x0 = loc[0] == 0 ? gd.x[0] : gd.xh[0];
        const TF y0 = loc[1] == 0 ? gd.y[0] : gd.yh[0];
        const std::vector<TF>& z = (loc[2] == 0) ? gd.z : gd.zh;
        const std::vector<TF>& dzi = (loc[2] == 0) ? gd.dzhi : gd.dzi;

        plk::calc_interpolation_factors_h(il.data(), fx.data(), xp.data(), x0, gd.dxi, n_particles);
        plk::calc_interpolation_factors_h(jl.data(), fy.data(), yp.data(), y0, gd.dyi, n_particles);
        plk::calc_interpolation_factors_v(kl.data(), fz.data(), zp.data(), z.data(), dzi.data(), loc[2], n_particles, gd.kcells);

        plk::diagnose_velocity(
            velocity.data(),
            fld_3d.data(),
            il.data(),
            jl.data(),
            kl.data(),
            fx.data(),
            fy.data(),
            fz.data(),
            n_particles,
            gd.jstride,
            gd.kstride);

        plk::add_tendency(
            tendency.data(),
            velocity.data(),
            n_particles);
    };

    diagnose_tendency(up, xpt, fields.mp.at("u")->fld, {1,0,0});
    diagnose_tendency(vp, ypt, fields.mp.at("v")->fld, {0,1,0});
    diagnose_tendency(wp, zpt, fields.mp.at("w")->fld, {0,0,1});
}


template<typename TF>
void Particle_lagrangian<TF>::integrate(Timeloop<TF>& timeloop)
{
    if (!sw_particle)
        return;

    auto& gd = grid.get_grid_data();

    // Integrate particle location with RK3/4 scheme.
    timeloop.exec(xp, xpt);
    timeloop.exec(yp, ypt);
    timeloop.exec(zp, zpt);

    // Quick hack: bounce particles from domain bottom/top.
    for (int n=0; n<n_particles; ++n)
        if (zp[n] < 0) zp[n] = -zp[n];

    // More quick hack: cyclic boundaries.
    for (int n=0; n<n_particles; ++n)
    {
        if (xp[n] >= gd.xsize)
            xp[n] -= gd.xsize;

        if (xp[n] < 0)
            xp[n] += gd.xsize;

        if (yp[n] >= gd.ysize)
            yp[n] -= gd.ysize;

        if (yp[n] < 0)
            yp[n] += gd.ysize;
    }
}
#endif


namespace
{
    template<typename TF>
    void read_coordinate(
        hid_t file_id, const char* dset_name,
        TF* buffer,
        hsize_t start,
        hsize_t count)
    {
        // Open dataset.
        hid_t dset = H5Dopen(file_id, dset_name, H5P_DEFAULT);
        hid_t filespace = H5Dget_space(dset);

        // Select hyperslab.
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);

        // Local memory layout is simple; continous 1D array of size `count`.
        hid_t memspace = H5Screate_simple(1, &count, NULL);

        // Setup collective IO where all tasks participate.
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // Read data.
        H5Dread(dset, get_hdf5_type<TF>(), memspace, filespace, plist_id, buffer);

        // Cleanup!
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Dclose(dset);
    }


    template <typename TF>
    void read_particles_parallel(
        const std::string& filename,
        std::vector<TF>& x,
        std::vector<TF>& y,
        std::vector<TF>& z,
        const int mpiid,
        const int nprocs)
    {
        // Open file with parallel HDF5.
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
        H5Pclose(plist_id);

        // Get total number of particles.
        hid_t dset_x = H5Dopen(file_id, "/x", H5P_DEFAULT);
        hid_t dspace = H5Dget_space(dset_x);
        hsize_t n_total;
        H5Sget_simple_extent_dims(dspace, &n_total, NULL);
        H5Sclose(dspace);

        // Each task has `ceil(n_total / nprocs)` particles (except last task, see next code block).
        const int np_per_task = std::ceil(TF(n_total) / nprocs);
        hsize_t start = mpiid * np_per_task;
        hsize_t count = np_per_task;

        // Last MPI task might have less if `n_total % nprocs != 0`.
        if (start + count > n_total)
            count = n_total - start;

        //std::cout << "n_total=" << n_total << ", n_per_task=" << np_per_task << ", mpiid=" << mpiid << ", start=" << start << ", count=" << count << std::endl;

        // Resize local vectors.
        // Use local size without buffer; these particles are only
        // temporarely on this task and shipped elsewhwere soon.
        x.resize(count);
        y.resize(count);
        z.resize(count);

        // Read each coordinate with hyperslab selection
        read_coordinate(file_id, "/x", x.data(), start, count);
        read_coordinate(file_id, "/y", y.data(), start, count);
        read_coordinate(file_id, "/z", z.data(), start, count);

        // Cleanup!
        H5Dclose(dset_x);
        H5Fclose(file_id);
    }
}


template<typename TF>
void Particle_lagrangian<TF>::load(const std::string& sim_name, const int iotime)
{
    if (!sw_particle)
        return;

    auto& md = master.get_MPI_data();

    std::ostringstream file_in;
    file_in << "particles." << std::setfill('0') << std::setw(7) << iotime << ".h5";

    #ifdef USEMPI
    read_particles_parallel<TF>(file_in.str(), xp, yp, zp, md.mpiid, md.nprocs);
    #else
    // TODO.
    // MPI tasks 0 reads and distributes data.
    //if (md.mpiid == 0)
    //{
    //    std::ostringstream file_in;
    //    file_in << "particles." << std::setfill('0') << std::setw(7) << iotime << ".h5";
    //    Hdf5_file h5_file_in(file_in.str(), Hdf5_mode::Read);

    //    Hdf5_variable<int> var_uid(h5_file_in, "particle_id");
    //    Hdf5_variable<TF> var_x(h5_file_in, "x");
    //    Hdf5_variable<TF> var_y(h5_file_in, "y");
    //    Hdf5_variable<TF> var_z(h5_file_in, "z");

    //    auto uid_in = var_uid.read();
    //    auto x_in = var_x.read();
    //    auto y_in = var_y.read();
    //    auto z_in = var_z.read();

    //    //file_in.close();


    //    // TEST TEST TEST: write back.
    //    std::ostringstream file_out;
    //    file_out << "particles_out." << std::setfill('0') << std::setw(7) << iotime << ".h5";
    //    Hdf5_file h5_file_out(file_out.str(), Hdf5_mode::Write);
    //    h5_file_out.add_dimension("particle_id", n_particles);

    //    Hdf5_variable<int> var_uid_out(h5_file_out, "particle_id", {"particle_id"});
    //    Hdf5_variable<TF> var_x_out(h5_file_out, "x", {"particle_id"});
    //    Hdf5_variable<TF> var_y_out(h5_file_out, "y", {"particle_id"});
    //    Hdf5_variable<TF> var_z_out(h5_file_out, "z", {"particle_id"});

    //    var_uid_out.insert(uid_in);
    //    var_x_out.insert(x_in);
    //    var_y_out.insert(y_in);
    //    var_z_out.insert(z_in);

    //    //file_out.close();

    //    throw 1;



    //    // No MPI; in data stays local.
    //    uid = uid_in;
    //    xp = x_in;
    //    yp = y_in;
    //    zp = z_in;

    //    up.resize(n_particles);
    //    vp.resize(n_particles);
    //    wp.resize(n_particles);

    //    xpt.resize(n_particles);
    //    ypt.resize(n_particles);
    //    zpt.resize(n_particles);

    //    il.resize(n_particles);
    //    jl.resize(n_particles);
    //    kl.resize(n_particles);

    //    fx.resize(n_particles);
    //    fy.resize(n_particles);
    //    fz.resize(n_particles);
    //}
    #endif

    throw 1;
}


template<typename TF>
void Particle_lagrangian<TF>::save(const std::string& sim_name, const int iotime)
{
    if (!sw_particle)
        return;

    // TODO, no restarts for now!
}


template<typename TF>
void Particle_lagrangian<TF>::create(Timeloop<TF>& timeloop)
{
    if (!sw_particle)
        return;
}


template<typename TF>
unsigned long Particle_lagrangian<TF>::get_time_limit(const unsigned long itime)
{
    if (!sw_dump)
        return Constants::ulhuge;

    return isampletime_dump - itime % isampletime_dump;
}


template<typename TF>
bool Particle_lagrangian<TF>::do_dump(const unsigned long itime)
{
    if (!sw_dump || itime % isampletime_dump != 0)
        return false;

    return true;
}


template<typename TF>
void Particle_lagrangian<TF>::dump(const int iotime)
{
    master.print_message("Saving raw particle dump\n");

    char file_name[256];
    std::sprintf(file_name, "particles.%07d", iotime);
    FILE* file = fopen(file_name, "wbx");

    // Check file opening and reading.
    bool success = (file != nullptr);

    if (success)
    {
        fwrite(xp.data(), sizeof(TF), n_particles, file);
        fwrite(yp.data(), sizeof(TF), n_particles, file);
        fwrite(zp.data(), sizeof(TF), n_particles, file);
    }

    if (!success)
    {
        #ifdef USEMPI
        std::cout << "SINGLE PROCESS EXCEPTION: saving particle dump " << file_name << " failed." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        #else
        throw std::runtime_error("ERROR: saving particle dump failed");
        #endif
    }

    fclose(file);
}


#ifdef FLOAT_SINGLE
template class Particle_lagrangian<float>;
#else
template class Particle_lagrangian<double>;
#endif
