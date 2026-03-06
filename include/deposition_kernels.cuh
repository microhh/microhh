/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef DEPOSITION_KERNELS_CUH
#define DEPOSITION_KERNELS_CUH

#include "constants.h"

namespace Deposition_kernels_g
{
    template<typename TF> __global__
    void calc_vd_water(
            TF* const __restrict__ fld,
            const TF* const __restrict__ ra,
            const TF* const __restrict__ ustar,
            const int* const __restrict__ water_mask,
            const TF diff_scl,
            const TF rwat,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*jstride;

            if (water_mask[ij] == 1)
            {
                const TF rb = TF(1) / (Constants::kappa<TF> * ustar[ij]) * diff_scl;
                fld[ij] = TF(1) / (ra[ij] + rb + rwat);
            }
        }
    }


    template<typename TF> __global__
    void calc_deposition_veg(
            TF* const __restrict__ vdo3,
            TF* const __restrict__ vdno,
            TF* const __restrict__ vdno2,
            TF* const __restrict__ vdhno3,
            TF* const __restrict__ vdh2o2,
            TF* const __restrict__ vdrooh,
            TF* const __restrict__ vdhcho,
            const TF* const __restrict__ lai,
            const TF* const __restrict__ rs,
            const TF* const __restrict__ ra,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ fraction,
            const TF* const __restrict__ rmes,
            const TF* const __restrict__ rsoil,
            const TF* const __restrict__ rcut,
            const TF* const __restrict__ diff_scl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*jstride;

            if (fraction[ij] < TF(1e-12))
                return;

            const TF hc = TF(10);
            const TF ra_inc = TF(14) * hc * lai[ij] / ustar[ij];
            const TF rb_fac = TF(2) / (Constants::kappa<TF> * ustar[ij]);

            // Rmes for NO and NO2 requires multiplication with rs (Ganzeveld et al. 1995).
            const TF rb_o3 = rb_fac * diff_scl[0];
            const TF rc_o3 = TF(1) / (TF(1) / (diff_scl[0] + rs[ij] + rmes[0]) + TF(1) / rcut[0] + TF(1) / (ra_inc + rsoil[0]));
            vdo3[ij] = TF(1) / (ra[ij] + rb_o3 + rc_o3);

            const TF rb_no = rb_fac * diff_scl[1];
            const TF rc_no = TF(1) / (TF(1) / (diff_scl[1] + rs[ij] + rmes[1]*rs[ij]) + TF(1) / rcut[1] + TF(1) / (ra_inc + rsoil[1]));
            vdno[ij] = TF(1) / (ra[ij] + rb_no + rc_no);

            const TF rb_no2 = rb_fac * diff_scl[2];
            const TF rc_no2 = TF(1) / (TF(1) / (diff_scl[2] + rs[ij] + rmes[2]*rs[ij]) + TF(1) / rcut[2] + TF(1) / (ra_inc + rsoil[2]));
            vdno2[ij] = TF(1) / (ra[ij] + rb_no2 + rc_no2);

            const TF rb_hno3 = rb_fac * diff_scl[3];
            const TF rc_hno3 = TF(1) / (TF(1) / (diff_scl[3] + rs[ij] + rmes[3]) + TF(1) / rcut[3] + TF(1) / (ra_inc + rsoil[3]));
            vdhno3[ij] = TF(1) / (ra[ij] + rb_hno3 + rc_hno3);

            const TF rb_h2o2 = rb_fac * diff_scl[4];
            const TF rc_h2o2 = TF(1) / (TF(1) / (diff_scl[4] + rs[ij] + rmes[4]) + TF(1) / rcut[4] + TF(1) / (ra_inc + rsoil[4]));
            vdh2o2[ij] = TF(1) / (ra[ij] + rb_h2o2 + rc_h2o2);

            const TF rb_rooh = rb_fac * diff_scl[5];
            const TF rc_rooh = TF(1) / (TF(1) / (diff_scl[5] + rs[ij] + rmes[5]) + TF(1) / rcut[5] + TF(1) / (ra_inc + rsoil[5]));
            vdrooh[ij] = TF(1) / (ra[ij] + rb_rooh + rc_rooh);

            const TF rb_hcho = rb_fac * diff_scl[6];
            const TF rc_hcho = TF(1) / (TF(1) / (diff_scl[6] + rs[ij] + rmes[6]) + TF(1) / rcut[6] + TF(1) / (ra_inc + rsoil[6]));
            vdhcho[ij] = TF(1) / (ra[ij] + rb_hcho + rc_hcho);
        }
    }


    template<typename TF> __global__
    void calc_deposition_soil(
            TF* const __restrict__ vdo3,
            TF* const __restrict__ vdno,
            TF* const __restrict__ vdno2,
            TF* const __restrict__ vdhno3,
            TF* const __restrict__ vdh2o2,
            TF* const __restrict__ vdrooh,
            TF* const __restrict__ vdhcho,
            const TF* const __restrict__ ra,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ fraction,
            const TF* const __restrict__ rsoil,
            const TF* const __restrict__ diff_scl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*jstride;

            if (fraction[ij] < TF(1e-12))
                return;

            const TF rb_fac = TF(1) / (Constants::kappa<TF> * ustar[ij]);

            vdo3[ij] = TF(1) / (ra[ij] + rb_fac * diff_scl[0] + rsoil[0]);
            vdno[ij] = TF(1) / (ra[ij] + rb_fac * diff_scl[1] + rsoil[1]);
            vdno2[ij] = TF(1) / (ra[ij] + rb_fac * diff_scl[2] + rsoil[2]);
            vdhno3[ij] = TF(1) / (ra[ij] + rb_fac * diff_scl[3] + rsoil[3]);
            vdh2o2[ij] = TF(1) / (ra[ij] + rb_fac * diff_scl[4] + rsoil[4]);
            vdrooh[ij] = TF(1) / (ra[ij] + rb_fac * diff_scl[5] + rsoil[5]);
            vdhcho[ij] = TF(1) / (ra[ij] + rb_fac * diff_scl[6] + rsoil[6]);
        }
    }


    template<typename TF> __global__
    void calc_deposition_wet(
            TF* const __restrict__ vdo3,
            TF* const __restrict__ vdno,
            TF* const __restrict__ vdno2,
            TF* const __restrict__ vdhno3,
            TF* const __restrict__ vdh2o2,
            TF* const __restrict__ vdrooh,
            TF* const __restrict__ vdhcho,
            const TF* const __restrict__ lai,
            const TF* const __restrict__ c_veg,
            const TF* const __restrict__ rs_wet,
            const TF* const __restrict__ rs_veg,
            const TF* const __restrict__ ra,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ fraction,
            const TF* const __restrict__ rmes,
            const TF* const __restrict__ rsoil,
            const TF* const __restrict__ rws,
            const TF* const __restrict__ diff_scl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*jstride;

            if (fraction[ij] < TF(1e-12))
                return;

            const TF hc = TF(10);
            const TF ra_inc = TF(14) * hc * lai[ij] / ustar[ij];
            const TF rb_fac = TF(1) / (Constants::kappa<TF> * ustar[ij]);
            const TF cveg = c_veg[ij];
            const TF one_min_cveg = TF(1) - cveg;

            // rcut replaced by rws for wet skin uptake. rb_veg == rb_soil, merged into rb.
            // Rmes for NO and NO2 requires multiplication with rs_wet (Ganzeveld et al. 1995).
            const TF rb_o3 = rb_fac * diff_scl[0];
            const TF rc_o3 = TF(1) / (TF(1) / (diff_scl[0] + rs_veg[ij] + rmes[0]) + TF(1) / rws[0] + TF(1) / (ra_inc + rsoil[0]));
            vdo3[ij] = cveg / (ra[ij] + rb_o3 + rc_o3) + one_min_cveg / (ra[ij] + rb_o3 + rsoil[0]);

            const TF rb_no = rb_fac * diff_scl[1];
            const TF rc_no = TF(1) / (TF(1) / (diff_scl[1] + rs_veg[ij] + rmes[1]*rs_wet[ij]) + TF(1) / rws[1] + TF(1) / (ra_inc + rsoil[1]));
            vdno[ij] = cveg / (ra[ij] + rb_no + rc_no) + one_min_cveg / (ra[ij] + rb_no + rsoil[1]);

            const TF rb_no2 = rb_fac * diff_scl[2];
            const TF rc_no2 = TF(1) / (TF(1) / (diff_scl[2] + rs_veg[ij] + rmes[2]*rs_wet[ij]) + TF(1) / rws[2] + TF(1) / (ra_inc + rsoil[2]));
            vdno2[ij] = cveg / (ra[ij] + rb_no2 + rc_no2) + one_min_cveg / (ra[ij] + rb_no2 + rsoil[2]);

            const TF rb_hno3 = rb_fac * diff_scl[3];
            const TF rc_hno3 = TF(1) / (TF(1) / (diff_scl[3] + rs_veg[ij] + rmes[3]) + TF(1) / rws[3] + TF(1) / (ra_inc + rsoil[3]));
            vdhno3[ij] = cveg / (ra[ij] + rb_hno3 + rc_hno3) + one_min_cveg / (ra[ij] + rb_hno3 + rsoil[3]);

            const TF rb_h2o2 = rb_fac * diff_scl[4];
            const TF rc_h2o2 = TF(1) / (TF(1) / (diff_scl[4] + rs_veg[ij] + rmes[4]) + TF(1) / rws[4] + TF(1) / (ra_inc + rsoil[4]));
            vdh2o2[ij] = cveg / (ra[ij] + rb_h2o2 + rc_h2o2) + one_min_cveg / (ra[ij] + rb_h2o2 + rsoil[4]);

            const TF rb_rooh = rb_fac * diff_scl[5];
            const TF rc_rooh = TF(1) / (TF(1) / (diff_scl[5] + rs_veg[ij] + rmes[5]) + TF(1) / rws[5] + TF(1) / (ra_inc + rsoil[5]));
            vdrooh[ij] = cveg / (ra[ij] + rb_rooh + rc_rooh) + one_min_cveg / (ra[ij] + rb_rooh + rsoil[5]);

            const TF rb_hcho = rb_fac * diff_scl[6];
            const TF rc_hcho = TF(1) / (TF(1) / (diff_scl[6] + rs_veg[ij] + rmes[6]) + TF(1) / rws[6] + TF(1) / (ra_inc + rsoil[6]));
            vdhcho[ij] = cveg / (ra[ij] + rb_hcho + rc_hcho) + one_min_cveg / (ra[ij] + rb_hcho + rsoil[6]);
        }
    }

} // namespace Deposition_kernels_g

#endif
