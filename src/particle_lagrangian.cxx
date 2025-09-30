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

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "constants.h"
#include "hdf5_interface.h"

#include "particle_lagrangian.h"

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


    template<typename TF>
    void calc_interpolation_factors_h(
        int* const restrict index,
        TF* const restrict factor,
        const TF* const restrict xp,
        const TF x0,
        const TF dxi,
        const int n_particles)
    {
        // For each particle, find index left of value, and calculate interpolation factor.
        // Equidistant grid in the horizontal, so index can be found directly.

        for (int i=0; i<n_particles; ++i)
        {
            const TF fi = (xp[i] - x0) * dxi;
            index[i] = static_cast<int>(fi);
            factor[i] = fi - index[i];
        }
    }


    template<typename TF>
    void calc_interpolation_factors_v(
        int* const restrict index,
        TF* const restrict factor,
        const TF* const restrict zp,
        const TF* const restrict z,
        const TF* const restrict dzi,
        const bool is_half_level,
        const int n_particles,
        const int kcells)
    {
        // For each particle, find index left of value, and calculate interpolation factor.
        // Non-equidistant grid in the vertical, so requires search using `std::upper_bound`.

        // This is slightly annoying...
        // For half levels, the spacing from zh[k] to zh[k+1] = dz[k]
        // For full levels, the spacing from z[k] to z[k+1] = dzh[k+1]
        const int dk = is_half_level ? 0 : 1;

        for (int i=0; i<n_particles; ++i)
        {
            // Use `upper_bound`; our `zh[0]` and `zh[1]` are both zero!
            const TF* it = std::upper_bound(z, z+kcells, zp[i]);
            const int k0 = static_cast<int>(it - z) - 1;

            index[i] = k0;
            factor[i] = (zp[i] - z[k0]) * dzi[k0+dk];
        }
    }


    template<typename TF>
    void diagnose_velocity(
        TF* const restrict vel_p,
        const TF* const restrict vel_3d,
        const int* const restrict il,
        const int* const restrict jl,
        const int* const restrict kl,
        const TF* const restrict fx,
        const TF* const restrict fy,
        const TF* const restrict fz,
        const int n_particles,
        const int jstride,
        const int kstride)
    {
        // Diagnose particle velocity by tri-linear interpolation of Eulerian velocity field to particle location.
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int n=0; n<n_particles; ++n)
        {
            const int ijk = il[n] + jl[n]*jstride + kl[n]*kstride;

            const TF fx1 = fx[n];
            const TF fy1 = fy[n];
            const TF fz1 = fz[n];

            const TF fx0 = TF(1) - fx1;
            const TF fy0 = TF(1) - fy1;
            const TF fz0 = TF(1) - fz1;

            vel_p[n] =
                fx0 * fy0 * fz0 * vel_3d[ijk               ] +
                fx1 * fy0 * fz0 * vel_3d[ijk + ii          ] +
                fx0 * fy1 * fz0 * vel_3d[ijk + jj          ] +
                fx0 * fy0 * fz1 * vel_3d[ijk + kk          ] +
                fx1 * fy1 * fz0 * vel_3d[ijk + ii + jj     ] +
                fx1 * fy0 * fz1 * vel_3d[ijk + ii + kk     ] +
                fx0 * fy1 * fz1 * vel_3d[ijk + jj + kk     ] +
                fx1 * fy1 * fz1 * vel_3d[ijk + ii + jj + kk];
        }
    }


    template<typename TF>
    void add_tendency(
        TF* const restrict tend_p,
        const TF* const restrict vel_p,
        const int n_particles)
    {
        // Add velocity to tendency.
        for (int n=0; n<n_particles; ++n)
            tend_p[n] += vel_p[n];
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

        calc_interpolation_factors_h(il.data(), fx.data(), xp.data(), x0, gd.dxi, n_particles);
        calc_interpolation_factors_h(jl.data(), fy.data(), yp.data(), y0, gd.dyi, n_particles);
        calc_interpolation_factors_v(kl.data(), fz.data(), zp.data(), z.data(), dzi.data(), loc[2], n_particles, gd.kcells);

        diagnose_velocity(
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

        add_tendency(
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


template<typename TF>
void Particle_lagrangian<TF>::load(const std::string& sim_name, const int iotime)
{
    if (!sw_particle)
        return;

    auto& md = master.get_MPI_data();

    // MPI tasks 0 reads and distributes data.
    if (md.mpiid == 0)
    {
        std::ostringstream oss;
        oss << "particles." << std::setfill('0') << std::setw(7) << iotime << ".h5";
        std::string file_name = oss.str();
        Hdf5_file h5_file(file_name, Hdf5_mode::Read);

        Hdf5_variable<int> var_uid(h5_file, "uid");
        Hdf5_variable<TF> var_x(h5_file, "x");
        Hdf5_variable<TF> var_y(h5_file, "y");
        Hdf5_variable<TF> var_z(h5_file, "z");

        auto uid_in = var_uid.read();
        auto x_in = var_x.read();
        auto y_in = var_y.read();
        auto z_in = var_z.read();

        // No MPI; in data stays local.
        uid = uid_in;
        xp = x_in;
        yp = y_in;
        zp = z_in;

        up.resize(n_particles);
        vp.resize(n_particles);
        wp.resize(n_particles);

        xpt.resize(n_particles);
        ypt.resize(n_particles);
        zpt.resize(n_particles);

        il.resize(n_particles);
        jl.resize(n_particles);
        kl.resize(n_particles);

        fx.resize(n_particles);
        fy.resize(n_particles);
        fz.resize(n_particles);
    }
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
