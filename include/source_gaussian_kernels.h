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

#ifndef SOURCE_GAUSSIAN_KERNELS_H
#define SOURCE_GAUSSIAN_KERNELS_H

#include "fast_math.h"
#include "constants.h"

namespace fm = Fast_math;

namespace Source_gaussian_kernels
{
    template<typename TF>
    std::vector<int> calc_shape(
            const TF* restrict x, const TF x0, const TF sigma_x, const TF line_x, int istart, int iend)
    {
        std::vector<int> range(2);
    
        int i = istart;
        range[0] = iend;
    
        for (; i<iend; ++i)
        {
            if ( x[i]-x0 + TF(4)*sigma_x > TF(0) )
            {
                range[0] = i;
                break;
            }
        }
    
        i = istart;
        for (; i<iend; ++i)
        {
            range[1] = iend;
    
            if ( x[i]-x0-line_x - TF(4)*sigma_x > TF(0) )
            {
                range[1] = i;
                break;
            }
        }
    
        return range;
    }


    template<typename TF>
    std::vector<int> calc_shape_profile(
            const TF* restrict emission_profile, const TF* restrict z, const int kstart, const int kend)
    {
        std::vector<int> range(2);

        for (int k=kstart; k<kend; ++k)
            if (emission_profile[k] > Constants::dsmall)
            {
                range[0] = k;
                break;
            }

        for (int k=kend-1; k>=kstart; --k)
            if (emission_profile[k] > Constants::dsmall)
            {
                range[1] = k+1;
                break;
            }

        return range;
    }


    template<typename TF, bool use_profile>
    void calc_source(
            TF* const __restrict__ st,
            const TF* const restrict x,
            const TF* const restrict y,
            const TF* const restrict z,
            const TF* const restrict emission_profile,
            const TF x0, const TF sigma_x, const TF line_x,
            const TF y0, const TF sigma_y, const TF line_y,
            const TF z0, const TF sigma_z, const TF line_z,
            const TF strength, const TF norm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        if (iend == istart || jend == jstart || kend == kstart)
            return;

        for(int k = kstart; k<kend; ++k)
            for(int j = jstart; j<jend; ++j)
                #pragma ivdep
                for(int i = istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    if (line_x != 0)
                    {
                        if (x[i] >= x0+line_x)
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else if (x[i] <= x0)
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0)       /fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    }
                    else if (line_y != 0)
                    {
                        if (y[j] >= y0+line_y)
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else if (y[j] <= y0)
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0)       /fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    }
                    else if (line_z != 0)
                    {
                        if (z[k] >= z0+line_z)
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else if (z[k] <= z0)
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0)       /fm::pow2(sigma_z));
                        else
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y));
                    }
                    else
                    {
                        if (use_profile)
                            st[ijk] += strength/norm*(exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y))
                                    * emission_profile[k]);
                        else
                            st[ijk] += strength/norm*exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    }
                }
        }
}
#endif
