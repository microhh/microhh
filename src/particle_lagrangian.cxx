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
#include "stats.h"

#include "particle_lagrangian.h"
#include "particle_lagrangian_kernels.h"
#include "particle_lagrangian_io.h"

namespace plk = Particle_lagrangian_kernels;
namespace plio = Particle_lagrangian_io;


template<typename TF>
Particle_lagrangian<TF>::Particle_lagrangian(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_particle = inputin.get_item<bool>("particle_lagrangian", "sw_particle", "", false);

    if (sw_particle)
    {
        // Particle properties. New properties only need to be added here; restart files,
        // MPI exchange, resizing, and time integration are handled by `Particle_data`.
        particles.add("x", Particle_kind::Prognostic, "Particle x-location", "m");
        particles.add("y", Particle_kind::Prognostic, "Particle y-location", "m");
        particles.add("z", Particle_kind::Prognostic, "Particle z-location", "m");

        particles.add("u", Particle_kind::Diagnostic, "Particle u-velocity", "m s-1");
        particles.add("v", Particle_kind::Diagnostic, "Particle v-velocity", "m s-1");
        particles.add("w", Particle_kind::Diagnostic, "Particle w-velocity", "m s-1");

        // Raw dump of all particles.
        sw_dump = inputin.get_item<bool>("particle_lagrangian", "sw_dump", "", false);

        if (sw_dump)
        {
            const int sampletime = inputin.get_item<int>("particle_lagrangian", "sampletime_dump", "");
            isampletime_dump = convert_to_itime(sampletime);
            szip_compression = inputin.get_item<int>("particle_lagrangian", "szip_compression", "", 0);
            dumplist = inputin.get_list<std::string>("particle_lagrangian", "dumplist", "");

            if (dumplist.empty())
                throw std::runtime_error("Particle dump enabled, but dumplist is empty");
            dump_file_id = -1;
        }

        // Profile and time serie statistics.
        sw_stats = inputin.get_item<bool>("particle_lagrangian", "sw_stats", "", false);

        // "Oversize" ratio arrays.
        reserve_ratio = inputin.get_item<TF>("particle_lagrangian", "reserve_ratio", "", 1.25);

        // Delayed start time of particles.
        const int starttime = inputin.get_item<int>("particle_lagrangian", "starttime", "", 0);
        istarttime = convert_to_itime(starttime);
    }
}


template<typename TF>
Particle_lagrangian<TF>::~Particle_lagrangian()
{
    if (sw_particle && sw_dump && dump_file_id >= 0)
        H5Fclose(dump_file_id);
}


template<typename TF>
void Particle_lagrangian<TF>::create(Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    if (!sw_particle)
        return;

    if (sw_dump)
    {
        const int iotime = timeloop.get_iotime();
        std::ostringstream file_out;
        file_out << "particle_dump." << std::setfill('0') << std::setw(7) << iotime << ".h5";

        for (const std::string& name : dumplist)
            if (!particles.has_property(name))
                throw std::runtime_error("Particle property \"" + name + "\" in dumplist does not exist");

        // Create single HDF5 file for particle dumps.
        plio::create_particle_dump<TF>(
            file_out.str(),
            dump_file_id,
            n_particles,
            szip_compression,
            dumplist,
            particles,
            master);
    }

    if (sw_stats && stats.get_switch())
    {
        std::string group_name = "particles";
        stats.add_time_series("mem_inc", "Count of increase in memory buffer size", "-", group_name);
        stats.add_time_series("mem_dec", "Count of decrease in memory buffer size", "-", group_name);
    }
}


#ifndef USECUDA
template<typename TF>
void Particle_lagrangian<TF>::exec(Timeloop<TF>& timeloop)
{
    // Calculate particle tendencies by tri-linear interpolation of Eulerian velocity fields to particle locations.
    if (!sw_particle || timeloop.get_itime() < istarttime)
        return;

    auto& gd = grid.get_grid_data();

    auto& xp = particles.get("x");
    auto& yp = particles.get("y");
    auto& zp = particles.get("z");

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

    diagnose_tendency(particles.get("u"), particles.get_tend("x"), fields.mp.at("u")->fld, {1,0,0});
    diagnose_tendency(particles.get("v"), particles.get_tend("y"), fields.mp.at("v")->fld, {0,1,0});
    diagnose_tendency(particles.get("w"), particles.get_tend("z"), fields.mp.at("w")->fld, {0,0,1});
}


template<typename TF>
void Particle_lagrangian<TF>::integrate(Timeloop<TF>& timeloop)
{
    if (!sw_particle || timeloop.get_itime() < istarttime)
        return;

    auto& gd = grid.get_grid_data();

    // Integrate all prognostic particle properties with RK3/4 scheme.
    for (const std::string& name : particles.get_names({Particle_kind::Prognostic}))
        timeloop.exec(particles.get(name), particles.get_tend(name));

    auto& xp = particles.get("x");
    auto& yp = particles.get("y");
    auto& zp = particles.get("z");

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
    int mem_inc_loc = 0;
    int mem_dec_loc = 0;

    plk::particle_exchange_parallel(
        particles,
        gd.xsize,
        gd.ysize,
        reserve_ratio,
        mem_inc_loc,
        mem_dec_loc,
        master);

    // Sum memory increase/decreases over all MPI ranks.
    master.sum(&mem_inc_loc, 1);
    master.sum(&mem_dec_loc, 1);

    mem_inc += mem_inc_loc;
    mem_dec += mem_dec_loc;

    // Resize interpolation indexes and factors.
    const int new_size = particles.size();

    adaptive_resize(il, new_size, reserve_ratio);
    adaptive_resize(jl, new_size, reserve_ratio);
    adaptive_resize(kl, new_size, reserve_ratio);

    adaptive_resize(fx, new_size, reserve_ratio);
    adaptive_resize(fy, new_size, reserve_ratio);
    adaptive_resize(fz, new_size, reserve_ratio);

    #else
    plk::particle_exchange_serial(
        xp,
        yp,
        gd.xsize,
        gd.ysize);
    #endif
}
#endif


template<typename TF>
void Particle_lagrangian<TF>::exec_stats(Stats<TF>& stats)
{
    if (!sw_particle || !sw_stats)
        return;

    stats.set_time_series("mem_inc", TF(mem_inc));
    stats.set_time_series("mem_dec", TF(mem_dec));
}


template<typename TF>
void Particle_lagrangian<TF>::load(const int iotime)
{
    if (!sw_particle)
        return;

    auto& md = master.get_MPI_data();
    auto& gd = grid.get_grid_data();

    std::ostringstream file_in;
    file_in << "particles." << std::setfill('0') << std::setw(7) << iotime << ".h5";

    // Read all prognostic and constant properties.
    const std::vector<std::string> names = particles.get_names({Particle_kind::Prognostic, Particle_kind::Constant});

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
    std::vector<std::vector<TF>> fields_in(names.size());

    // Read particles to their "home" task using parallel HDF5.
    plio::read_particles<TF>(file_in.str(), names, uid_in, get_pointers(fields_in), n_particles, master);

    // Send particles from "home" task to actual location in domain.
    auto get_index = [&](const std::string& name)
    {
        return std::find(names.begin(), names.end(), name) - names.begin();
    };

    plio::distribute_particles<TF>(
        particles.uid, particles.get_fields(names),
        uid_in, get_pointers(fields_in),
        fields_in[get_index("x")], fields_in[get_index("y")],
        gd.xsize, gd.ysize,
        master);

    #else
    /*
     * Read particles using serial HDF5.
     */
    plio::read_particles<TF>(file_in.str(), names, particles.uid, particles.get_fields(names), n_particles, master);
    #endif

    // Resize other properties and tendencies.
    const int n_local = particles.size();
    particles.resize(n_local);

    // Interpolation indexes and factors.
    il.resize(n_local);
    jl.resize(n_local);
    kl.resize(n_local);

    fx.resize(n_local);
    fy.resize(n_local);
    fz.resize(n_local);
}


template<typename TF>
void Particle_lagrangian<TF>::save(const int iotime)
{
    if (!sw_particle || iotime == 0)
        return;

    std::ostringstream file_out;
    file_out << "particles." << std::setfill('0') << std::setw(7) << iotime << ".h5";

    master.print_message("Saving \"%s\" ... ", file_out.str().c_str());

    // Save all prognostic and constant properties.
    const std::vector<std::string> names = particles.get_names({Particle_kind::Prognostic, Particle_kind::Constant});

    // Create new HDF5 file for each restart file.
    hid_t restart_file_id;
    plio::create_particle_restart<TF>(
        file_out.str(),
        restart_file_id,
        n_particles,
        names,
        particles,
        master);

    // Gather particles back to their original MPI tasks based on `uid`.
    std::vector<int> uid_local;
    std::vector<std::vector<TF>> fields_local(names.size());

    plio::gather_particles<TF>(
        uid_local, get_pointers(fields_local),
        particles.uid, particles.get_fields(names),
        n_particles,
        master);

    // Write gathered particles to restart file.
    auto& md = master.get_MPI_data();
    plio::write_particles_restart<TF>(
        restart_file_id,
        names,
        get_pointers(fields_local),
        n_particles,
        md.mpiid,
        md.nprocs);

    // Close the restart file.
    plio::h5_check(H5Fclose(restart_file_id), "closing " + file_out.str());

    master.print_message("OK\n");
}


template<typename TF>
unsigned long Particle_lagrangian<TF>::get_time_limit(const unsigned long itime)
{
    if (!sw_particle || !sw_dump)
        return Constants::ulhuge;

    return isampletime_dump - itime % isampletime_dump;
}


template<typename TF>
bool Particle_lagrangian<TF>::do_dump(const unsigned long itime)
{
    if (!sw_particle || !sw_dump || itime % isampletime_dump != 0)
        return false;

    return true;
}


template<typename TF>
void Particle_lagrangian<TF>::dump(Timeloop<TF>& timeloop)
{
    if (timeloop.get_itime() < istarttime)
        return;

    master.print_message("Saving particle dump for time %f\n", timeloop.get_time());

    auto& md = master.get_MPI_data();

    const int iotime = timeloop.get_iotime();
    const double time = timeloop.get_time();

    // Gather particles back to their original MPI tasks based on `uid`.
    std::vector<int> uid_local;
    std::vector<std::vector<TF>> fields_local(dumplist.size());

    plio::gather_particles<TF>(
        uid_local, get_pointers(fields_local),
        particles.uid, particles.get_fields(dumplist),
        n_particles,
        master);

    // Write particles to HDF5 file with parallel IO.
    plio::write_particles_dump<TF>(
        dump_file_id,
        dumplist,
        get_pointers(fields_local),
        time,
        n_particles,
        md.mpiid,
        md.nprocs);
}


#ifdef FLOAT_SINGLE
template class Particle_lagrangian<float>;
#else
template class Particle_lagrangian<double>;
#endif