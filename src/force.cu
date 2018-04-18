/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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

#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
#include "finite_difference.h"
#include "constants.h"
#include "tools.h"
#include "boundary.h"
#include "timeloop.h"
#include "data_block.h"
#include "force.h"

using namespace Finite_difference::O2;

namespace
{
    template<typename TF> __global__
    void flux_step_1_g(TF* const __restrict__ aSum, const TF* const __restrict__ a,
                       const TF* const __restrict__ dz,
                       const int jj, const int kk,
                       const int istart, const int jstart, const int kstart,
                       const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            aSum [ijk] = a[ijk]*dz[k];
        }
    }

    template<typename TF> __global__
    void flux_step_2_g(TF* const __restrict__ ut,
                       const TF fbody,
                       const int jj, const int kk,
                       const int istart, const int jstart, const int kstart,
                       const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += fbody;
        }
    }

    template<typename TF> __global__
    void coriolis_2nd_g(TF* const __restrict__ ut, TF* const __restrict__ vt,
                        TF* const __restrict__ u,  TF* const __restrict__ v,
                        TF* const __restrict__ ug, TF* const __restrict__ vg,
                        const TF fc, const TF ugrid, const TF vgrid,
                        const int jj, const int kk,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += fc * (0.25*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[k]);
            vt[ijk] -= fc * (0.25*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[k]);
        }
    }

    template<typename TF> __global__
    void coriolis_4th_g(TF* const __restrict__ ut, TF* const __restrict__ vt,
                        TF* const __restrict__ u,  TF* const __restrict__ v,
                        TF* const __restrict__ ug, TF* const __restrict__ vg,
                        const TF fc, const TF ugrid, const TF vgrid,
                        const int jj, const int kk,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend)
    {
        using namespace Finite_difference::O4;

        const int i   = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j   = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k   = blockIdx.z + kstart;
        const int ii  = 1;
        const int ii2 = 2;
        const int jj2 = 2*jj;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += fc * ( ( ci0*(ci0*v[ijk-ii2-jj ] + ci1*v[ijk-ii-jj ] + ci2*v[ijk-jj    ] + ci3*v[ijk+ii-jj  ])
                              + ci1*(ci0*v[ijk-ii2    ] + ci1*v[ijk-ii    ] + ci2*v[ijk       ] + ci3*v[ijk+ii     ])
                              + ci2*(ci0*v[ijk-ii2+jj ] + ci1*v[ijk-ii+jj ] + ci2*v[ijk+jj    ] + ci3*v[ijk+ii+jj  ])
                              + ci3*(ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii+jj2] + ci2*v[ijk+jj2   ] + ci3*v[ijk+ii+jj2 ]) )
                       + vgrid - vg[k] );

            vt[ijk] -= fc * ( ( ci0*(ci0*u[ijk-ii-jj2 ] + ci1*u[ijk-jj2   ] + ci2*u[ijk+ii-jj2] + ci3*u[ijk+ii2-jj2])
                              + ci1*(ci0*u[ijk-ii-jj  ] + ci1*u[ijk-jj    ] + ci2*u[ijk+ii-jj ] + ci3*u[ijk+ii2-jj ])
                              + ci2*(ci0*u[ijk-ii     ] + ci1*u[ijk       ] + ci2*u[ijk+ii    ] + ci3*u[ijk+ii2    ])
                              + ci3*(ci0*u[ijk-ii+jj  ] + ci1*u[ijk+jj    ] + ci2*u[ijk+ii+jj ] + ci3*u[ijk+ii2+jj ]) )
                       + ugrid - ug[k]);
        }
    }

    template<typename TF> __global__
    void advec_wls_2nd_g(TF* const __restrict__ st, TF* const __restrict__ s,
                         const TF* const __restrict__ wls, const TF* const __restrict__ dzhi,
                         const int istart, const int jstart, const int kstart,
                         const int iend,   const int jend,   const int kend,
                         const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            if (wls[k] > 0.)
                st[ijk] -=  wls[k] * (s[k]-s[k-1])*dzhi[k];
            else
                st[ijk] -=  wls[k] * (s[k+1]-s[k])*dzhi[k+1];
        }
    }

    template<typename TF> __global__
    void large_scale_source_g(TF* const __restrict__ st, TF* const __restrict__ sls,
                              const int istart, const int jstart, const int kstart,
                              const int iend,   const int jend,   const int kend,
                              const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            st[ijk] += sls[k];
        }
    }

    template<typename TF> __global__
    void nudging_tendency_g(TF* const __restrict__ st, TF* const __restrict__ smn,
			                TF* const __restrict__ snudge, TF* const __restrict__ nudge_fac,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend,
                            const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            st[ijk] += - nudge_fac[k] * (smn[k]-snudge[k]);

        }
    }

    template<typename TF> __global__
    void calc_time_dependent_prof_g(TF* const __restrict__ prof, const TF* const __restrict__ data,
                                    const double fac0, const double fac1,
                                    const int index0,  const int index1,
                                    const int kmax,    const int kgc)
    {
        const int k = blockIdx.x*blockDim.x + threadIdx.x;
        const int kk = kmax;

        if (k < kmax)
            prof[k+kgc] = fac0*data[index0*kk+k] + fac1*data[index1*kk+k];
    }
} // end namespace

template<typename TF>
void Force<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int nmemsize  = gd.kcells*sizeof(TF);

    if (swlspres== Large_scale_pressure_type::geo_wind)
    {
        cuda_safe_call(cudaMalloc(&ug_g, nmemsize));
        cuda_safe_call(cudaMalloc(&vg_g, nmemsize));

        cuda_safe_call(cudaMemcpy(ug_g, ug.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(vg_g, vg.data(), nmemsize, cudaMemcpyHostToDevice));

        if (tdep_geo.sw)
        {

            for (auto& it : tdep_geo.data_g)
            {
                int nmemsize2 = gd.kmax*tdep_geo.time[it.first].size()*sizeof(TF);
                cuda_safe_call(cudaMalloc(&it.second, nmemsize2));
                cuda_safe_call(cudaMemcpy(it.second, tdep_geo.data[it.first].data(), nmemsize2, cudaMemcpyHostToDevice));
            }
        }
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        for (auto& it : lslist)
        {
            cuda_safe_call(cudaMalloc(&lsprofs_g[it], nmemsize));
            cuda_safe_call(cudaMemcpy(lsprofs_g[it], lsprofs[it].data(), nmemsize, cudaMemcpyHostToDevice));
        }
        if (tdep_ls.sw)
        {

            for (auto& it : tdep_ls.data_g)
            {
                int nmemsize2 = gd.kmax*tdep_ls.time[it.first].size()*sizeof(TF);
                cuda_safe_call(cudaMalloc(&it.second, nmemsize2));
                cuda_safe_call(cudaMemcpy(it.second, tdep_ls.data[it.first].data(), nmemsize2, cudaMemcpyHostToDevice));
            }
        }
    }

    if (swnudge == Nudging_type::enabled)
    {
        for (auto& it : nudgelist)
        {
            cuda_safe_call(cudaMalloc(&nudgeprofs_g[it], nmemsize));
            cuda_safe_call(cudaMemcpy(nudgeprofs_g[it], nudgeprofs[it].data(), nmemsize, cudaMemcpyHostToDevice));
        }
        cuda_safe_call(cudaMalloc(&nudge_factor_g, nmemsize));
        cuda_safe_call(cudaMemcpy(nudge_factor_g, nudge_factor.data(), nmemsize, cudaMemcpyHostToDevice));

        if (tdep_nudge.sw)
        {

            for (auto& it : tdep_nudge.data_g)
            {
                int nmemsize2 = gd.kmax*tdep_nudge.time[it.first].size()*sizeof(TF);
                cuda_safe_call(cudaMalloc(&it.second, nmemsize2));
                cuda_safe_call(cudaMemcpy(it.second, tdep_nudge.data[it.first].data(), nmemsize2, cudaMemcpyHostToDevice));
            }
        }
    }

    if (swwls == Large_scale_subsidence_type::enabled)
    {
        cuda_safe_call(cudaMalloc(&wls_g, nmemsize));
        cuda_safe_call(cudaMemcpy(wls_g, wls.data(), nmemsize, cudaMemcpyHostToDevice));

        if (tdep_wls.sw)
        {
            int nmemsize2 = gd.kmax*tdep_geo.time["wls"].size()*sizeof(TF);
            cuda_safe_call(cudaMalloc(&tdep_geo.data_g["wls"], nmemsize2));
            cuda_safe_call(cudaMemcpy(tdep_geo.data_g["wls"], tdep_geo.data["wls"].data(), nmemsize2, cudaMemcpyHostToDevice));
        }
    }
}

template<typename TF>
void Force<TF>::clear_device()
{
    if (swlspres== Large_scale_pressure_type::geo_wind)
    {
        cuda_safe_call(cudaFree(ug_g));
        cuda_safe_call(cudaFree(vg_g));
        if (tdep_geo.sw)
        {
            for (auto& it : tdep_geo.data_g)
                cuda_safe_call(cudaFree(it.second));
        }
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        for(auto& it : lsprofs_g)
            cuda_safe_call(cudaFree(it.second));
            if (tdep_ls.sw)
            {
                for (auto& it : tdep_ls.data_g)
                    cuda_safe_call(cudaFree(it.second));
            }
    }

    if (swnudge == Nudging_type::enabled)
    {
        for(auto& it : nudgeprofs_g)
            cuda_safe_call(cudaFree(it.second));
            cuda_safe_call(cudaFree(nudge_factor_g));
            if (tdep_nudge.sw)
            {
                for (auto& it : tdep_nudge.data_g)
                    cuda_safe_call(cudaFree(it.second));
            }
    }

    if (swwls == Large_scale_subsidence_type::enabled)
    {
        cuda_safe_call(cudaFree(wls_g));
        if (tdep_wls.sw)
        {
            for (auto& it : tdep_wls.data_g)
                cuda_safe_call(cudaFree(it.second));
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Force<TF>::exec(double dt)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    if (swlspres == Large_scale_pressure_type::fixed_flux)
    {
        auto tmp = fields.get_tmp_g();

        TF uavg  = field3d_operators.calc_mean(fields.mp.at("u")->fld_g);
        TF utavg = field3d_operators.calc_mean(fields.mt.at("u")->fld_g);

        fields.release_tmp_g(tmp);

        const TF fbody = (uflux - uavg - grid.utrans) / dt - utavg;

        flux_step_2_g<TF><<<gridGPU, blockGPU>>>(
            fields.mt.at("u")->fld_g,
            fbody,
            gd.icells, gd.ijcells,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend,    gd.jend,   gd.kend);
        cuda_check_error();
    }
    else if (swlspres== Large_scale_pressure_type::geo_wind)
    {
        if (grid.swspatialorder == "2")
        {
            coriolis_2nd_g<<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
                ug_g, vg_g, fc, grid.utrans, grid.vtrans,
                gd.icells, gd.ijcells,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend);
            cuda_check_error();
        }
        else if (grid.swspatialorder == "4")
        {
            coriolis_4th_g<<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
                ug_g, vg_g, fc, grid.utrans, grid.vtrans,
                gd.icells, gd.ijcells,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend);
            cuda_check_error();
        }
    }

    if (swls == Large_scale_tendency_type::enabled)
    {
        for (auto& it : lslist)
        {
            large_scale_source_g<<<gridGPU, blockGPU>>>(
                fields.st.at(it)->fld_g, lsprofs_g.at(it),
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }

    if (swnudge == Nudging_type::enabled)
    {
        for (auto& it : nudgelist)
        {
            nudging_tendency_g<<<gridGPU, blockGPU>>>(
                fields.st.at(it)->fld_g, fields.sp.at(it)->fld_mean_g,
                nudgeprofs_g.at(it), nudge_factor_g,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }

    if (swwls == Large_scale_subsidence_type::enabled)
    {
        for (auto& it : fields.st)
        {
            advec_wls_2nd_g<<<gridGPU, blockGPU>>>(
                fields.st.at(it.first)->fld_g, fields.sp.at(it.first)->fld_mean_g, wls_g, gd.dzhi_g,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }
}
#endif

#ifdef USECUDA
template <typename TF>
void Force<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (tdep_ls.sw)
        update_time_dependent_profs_g(timeloop, lsprofs_g, tdep_ls);
    if (tdep_nudge.sw)
        update_time_dependent_profs_g(timeloop, nudgeprofs_g, tdep_nudge);
    if (tdep_geo.sw)
    {
        update_time_dependent_prof_g(timeloop, ug_g, tdep_geo, "u");
        update_time_dependent_prof_g(timeloop, vg_g, tdep_geo, "v");
    }
    if (tdep_wls.sw)
        update_time_dependent_prof_g(timeloop, wls_g, tdep_wls, "w");
}
#endif

#ifdef USECUDA
template<typename TF>
void Force<TF>::update_time_dependent_profs_g(Timeloop<TF>& timeloop, std::map<std::string, TF*> profiles, Time_dep timedep)
{
    auto& gd = grid.get_grid_data();
    const int blockk = 128;
    const int gridk  = gd.kmax/blockk + (gd.kmax%blockk > 0);

    // Loop over all profiles which might be time dependent
    for (auto& it : timedep.data_g)
    {
        // Get/calculate the interpolation indexes/factors
        int index0 = 0, index1 = 0;
        TF fac0 = 0., fac1 = 0.;

        timeloop.get_interpolation_factors(index0, index1, fac0, fac1, timedep.time[it.first]);

        // Calculate the new vertical profile
        calc_time_dependent_prof_g<<<gridk, blockk>>>(
            profiles[it.first], it.second, fac0, fac1, index0, index1, gd.kmax, gd.kgc);
        cuda_check_error();
    }
}
#endif

#ifdef USECUDA
template <typename TF>
void Force<TF>::update_time_dependent_prof_g(Timeloop<TF>& timeloop, TF* prof, Time_dep timedep, const std::string& name)
{
    auto& gd = grid.get_grid_data();
    const int blockk = 128;
    const int gridk  = gd.kmax/blockk + (gd.kmax%blockk > 0);

    // Get/calculate the interpolation indexes/factors
    int index0 = 0, index1 = 0;
    TF fac0 = 0., fac1 = 0.;
    timeloop.get_interpolation_factors(index0, index1, fac0, fac1, timedep.time[name]);

    // Calculate the new vertical profile
    calc_time_dependent_prof_g<<<gridk, blockk>>>(
        prof, timedep.data_g[name], fac0, fac1, index0, index1, gd.kmax, gd.kgc);
    cuda_check_error();
}
#endif

template class Force<double>;
template class Force<float>;
