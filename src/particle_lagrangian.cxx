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

#include "particle_lagrangian.h"
#include "particle_lagrangian_kernels.h"
#include "particle_lagrangian_io.h"

namespace plk = Particle_lagrangian_kernels;
namespace plio = Particle_lagrangian_io;

namespace
{

}


template<typename TF>
Particle_lagrangian<TF>::Particle_lagrangian(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_particle = inputin.get_item<bool>("particle_lagrangian", "sw_particle", "", false);

    if (sw_particle)
    {
        // Raw dump of all particles.
        sw_dump = inputin.get_item<bool>("particle_lagrangian", "sw_dump", "", false);

        if (sw_dump)
        {
            const int sampletime = inputin.get_item<int>("particle_lagrangian", "sampletime_dump", "");
            isampletime_dump = convert_to_itime(sampletime);
        }

        reserve_ratio = inputin.get_item<TF>("particle_lagrangian", "reserve_ratio", "", 1.25);
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

        plk::calc_interpolation_factors_h(il.data(), fx.data(), xp.data(), x0, gd.dxi, xp.size());
        plk::calc_interpolation_factors_h(jl.data(), fy.data(), yp.data(), y0, gd.dyi, xp.size());
        plk::calc_interpolation_factors_v(kl.data(), fz.data(), zp.data(), z.data(), dzi.data(), loc[2], xp.size(), gd.kcells);

        plk::diagnose_velocity(
            velocity.data(),
            fld_3d.data(),
            il.data(),
            jl.data(),
            kl.data(),
            fx.data(),
            fy.data(),
            fz.data(),
            xp.size(),
            gd.jstride,
            gd.kstride);

        plk::add_tendency(
            tendency.data(),
            velocity.data(),
            xp.size());
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
    // TODO: I think we need to manipulate tendencies to do this correctly...
    for (int n=0; n<xp.size(); ++n)
    {
        if (zp[n] < 0)
            zp[n] = -zp[n];

        if (zp[n] > gd.zsize)
            zp[n] = 2*gd.zsize - zp[n];
    }

    // Neighbour-neighbour exchange and cyclic boundary conditions.
    #ifdef USEMPI
    plk::particle_exchange_parallel(
        uid,
        xp,
        yp,
        zp,
        xpt,
        ypt,
        zpt,
        gd.xsize,
        gd.ysize,
        reserve_ratio,
        master);

    // Resize non-communicated vectors.
    const int new_size = xp.size();

    // These need to be included in the particle exchange
    // when the velocities become prognostic.
    plk::adaptive_resize(up, new_size, reserve_ratio);
    plk::adaptive_resize(vp, new_size, reserve_ratio);
    plk::adaptive_resize(wp, new_size, reserve_ratio);

    // These are always diagnostic.
    plk::adaptive_resize(il, new_size, reserve_ratio);
    plk::adaptive_resize(jl, new_size, reserve_ratio);
    plk::adaptive_resize(kl, new_size, reserve_ratio);

    plk::adaptive_resize(fx, new_size, reserve_ratio);
    plk::adaptive_resize(fy, new_size, reserve_ratio);
    plk::adaptive_resize(fz, new_size, reserve_ratio);

    #else
    plk::periodic_exchange_serial(
        xp,
        yp,
        gd.xsize,
        gd.ysize);
    #endif
}
#endif


template<typename TF>
void Particle_lagrangian<TF>::load(const std::string& sim_name, const int iotime)
{
    if (!sw_particle)
        return;

    auto& md = master.get_MPI_data();
    auto& gd = grid.get_grid_data();

    std::ostringstream file_in;
    file_in << "particles." << std::setfill('0') << std::setw(7) << iotime << ".h5";

    #ifdef USEMPI
    /*
     * Read particles in parallel using HDF5 and MPI.
     * Idea: assign each particle a "home" MPI task defined by:
     *     `mpiid = ceil(n_particles / nprocs)`
     * Since the particles are sorted in the input file, each MPI task
     * can read a continous block of data from `particles.0000000.h5` with
     * parallel HDF5, and then then send the particles to the
     * MPI tasks where they belong using `MPI_Alltoallv`.
     * This balances I/O across all MPI tasks for large amounts of particles.
     */

    // Read buffers.
    std::vector<int> uid_in;
    std::vector<TF> xp_in;
    std::vector<TF> yp_in;
    std::vector<TF> zp_in;

    // Read particles to their "home" task using parallel HDF5.
    plio::read_particles_parallel<TF>(file_in.str(), uid_in, xp_in, yp_in, zp_in, md.mpiid, md.nprocs);

    // Send particles from "home" task to actual location in domain.
    plio::distribute_particles(uid, xp, yp, zp, uid_in, xp_in, yp_in, zp_in,  gd.xsize, gd.ysize, master);

    #else
    /*
     * Read particles using simple serial HDF5.
     */
    plio::read_particles_serial(file_in.str(), uid, xp, yp, zp)
    #endif

    // Resize other properties.
    const int n_local = xp.size();

    // Velocities.
    up.resize(n_local);
    vp.resize(n_local);
    wp.resize(n_local);

    // Location tendencies.
    xpt.resize(n_local);
    ypt.resize(n_local);
    zpt.resize(n_local);

    // Interpolation indexes and factors.
    il.resize(n_local);
    jl.resize(n_local);
    kl.resize(n_local);

    fx.resize(n_local);
    fy.resize(n_local);
    fz.resize(n_local);
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
    //master.print_message("Saving raw particle dump\n");

    char file_name[256];
    std::sprintf(file_name, "particles.%07d", iotime);
    FILE* file = fopen(file_name, "wbx");

    // Check file opening and reading.
    bool success = (file != nullptr);

    if (success)
    {
        fwrite(xp.data(), sizeof(TF), xp.size(), file);
        fwrite(yp.data(), sizeof(TF), yp.size(), file);
        fwrite(zp.data(), sizeof(TF), zp.size(), file);
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