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

#ifndef PARTICLE_LAGRANGIAN_KERNELS_H
#define PARTICLE_LAGRANGIAN_KERNELS_H

namespace Particle_lagrangian_kernels
{



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
#endif
