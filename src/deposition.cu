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

#include "deposition.h"
#include "tools.h"
#include "constants.h"

namespace
{
    template<typename TF> __global__
    void calc_deposition_veg_g(
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
    void calc_deposition_soil_g(
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
    void calc_deposition_wet_g(
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


    template<typename TF> __global__
    void calc_vd_water_g(
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
    void calc_tiled_mean_g(
            TF* const __restrict__ fld_mean,
            const TF* const __restrict__ fld_veg,
            const TF* const __restrict__ fld_soil,
            const TF* const __restrict__ fld_wet,
            const TF* const __restrict__ tile_frac_veg,
            const TF* const __restrict__ tile_frac_soil,
            const TF* const __restrict__ tile_frac_wet,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;

            fld_mean[ij] = (
                tile_frac_veg [ij] * fld_veg [ij] +
                tile_frac_soil[ij] * fld_soil[ij] +
                tile_frac_wet [ij] * fld_wet [ij] );
        }
    }
}

#ifdef USECUDA
template <typename TF>
void Deposition<TF>::update_time_dependent(
        Timeloop<TF>& timeloop,
        Boundary<TF>& boundary,
        TF* restrict vdo3_g,
        TF* restrict vdno_g,
        TF* restrict vdno2_g,
        TF* restrict vdhno3_g,
        TF* restrict vdh2o2_g,
        TF* restrict vdrooh_g,
        TF* restrict vdhcho_g
        )
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // Get information from the land-surface model:
    auto& tiles = boundary.get_tiles();
    int* water_mask_g = boundary.get_water_mask_g();
    TF* lai_g = boundary.get_lai_g();
    TF* c_veg_g = boundary.get_c_veg_g();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 grid_gpu (gridi,  gridj,  1);
    dim3 block_gpu(blocki, blockj, 1);

    auto& dep_veg  = deposition_tiles.at("veg");
    auto& dep_soil = deposition_tiles.at("soil");
    auto& dep_wet  = deposition_tiles.at("wet");

    calc_deposition_veg_g<TF><<<grid_gpu, block_gpu>>>(
            dep_veg.vd_g.at("o3"),
            dep_veg.vd_g.at("no"),
            dep_veg.vd_g.at("no2"),
            dep_veg.vd_g.at("hno3"),
            dep_veg.vd_g.at("h2o2"),
            dep_veg.vd_g.at("rooh"),
            dep_veg.vd_g.at("hcho"),
            lai_g,
            tiles.at("veg").rs_g,
            tiles.at("veg").ra_g,
            tiles.at("veg").ustar_g,
            tiles.at("veg").fraction_g,
            rmes_g,
            rsoil_g,
            rcut_g,
            diff_scl_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    calc_deposition_soil_g<TF><<<grid_gpu, block_gpu>>>(
            dep_soil.vd_g.at("o3"),
            dep_soil.vd_g.at("no"),
            dep_soil.vd_g.at("no2"),
            dep_soil.vd_g.at("hno3"),
            dep_soil.vd_g.at("h2o2"),
            dep_soil.vd_g.at("rooh"),
            dep_soil.vd_g.at("hcho"),
            tiles.at("soil").ra_g,
            tiles.at("soil").ustar_g,
            tiles.at("soil").fraction_g,
            rsoil_g,
            diff_scl_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    calc_deposition_wet_g<TF><<<grid_gpu, block_gpu>>>(
            dep_wet.vd_g.at("o3"),
            dep_wet.vd_g.at("no"),
            dep_wet.vd_g.at("no2"),
            dep_wet.vd_g.at("hno3"),
            dep_wet.vd_g.at("h2o2"),
            dep_wet.vd_g.at("rooh"),
            dep_wet.vd_g.at("hcho"),
            lai_g,
            c_veg_g,
            tiles.at("wet").rs_g,
            tiles.at("veg").rs_g,
            tiles.at("wet").ra_g,
            tiles.at("wet").ustar_g,
            tiles.at("wet").fraction_g,
            rmes_g,
            rsoil_g,
            rws_g,
            diff_scl_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    auto calc_vd_g = [&](TF* vd, const std::string& sp)
    {
        calc_tiled_mean_g<TF><<<grid_gpu, block_gpu>>>(
                vd,
                dep_veg.vd_g.at(sp),
                dep_soil.vd_g.at(sp),
                dep_wet.vd_g.at(sp),
                tiles.at("veg").fraction_g,
                tiles.at("soil").fraction_g,
                tiles.at("wet").fraction_g,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        cuda_check_error();

        // Use wet-tile u* and ra: calculated in lsm with f_wet = 100%.
        const int s = species_idx.at(sp);
        calc_vd_water_g<TF><<<grid_gpu, block_gpu>>>(
                vd,
                tiles.at("wet").ra_g,
                tiles.at("wet").ustar_g,
                water_mask_g,
                diff_scl[s], rwat[s],
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        cuda_check_error();

        // TODO: spatial_avg_vd.
    };

    calc_vd_g(vdo3_g,  "o3");
    calc_vd_g(vdno_g,  "no");
    calc_vd_g(vdno2_g, "no2");
    calc_vd_g(vdhno3_g,"hno3");
    calc_vd_g(vdh2o2_g,"h2o2");
    calc_vd_g(vdrooh_g,"rooh");
    calc_vd_g(vdhcho_g,"hcho");
}

template <typename TF>
void Deposition<TF>::prepare_device()
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // Allocate GPU arrays.
    for (auto& tile : deposition_tiles)
        for (auto& [sp, idx] : species_idx)
            tile.second.vd_g[sp].allocate(gd.ijcells);

    const int size = rmes.size();
    rmes_g.allocate(size);
    rsoil_g.allocate(size);
    rcut_g.allocate(size);
    rws_g.allocate(size);
    rwat_g.allocate(size);
    diff_g.allocate(size);
    diff_scl_g.allocate(size);
    henry_g.allocate(size);
    f0_g.allocate(size);

    // Copy data to device.
    const int memsize_ij  = gd.ijcells * sizeof(TF);
    const int memsize_res = size * sizeof(TF);

    for (auto& tile : deposition_tiles)
        for (auto& [sp, idx] : species_idx)
            cuda_safe_call(cudaMemcpy(tile.second.vd_g.at(sp), tile.second.vd.at(sp).data(), memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(rmes_g,    rmes.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rsoil_g,   rsoil.data(),   memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rcut_g,    rcut.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rws_g,     rws.data(),     memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rwat_g,    rwat.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diff_g,    diff.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diff_scl_g,diff_scl.data(),memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(henry_g,   henry.data(),   memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(f0_g,      f0.data(),      memsize_res, cudaMemcpyHostToDevice));
}
#endif

#ifdef FLOAT_SINGLE
template class Deposition<float>;
#else
template class Deposition<double>;
#endif