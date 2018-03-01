/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "force.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
#include "finite_difference.h"
#include "constants.h"
#include "tools.h"
#include "boundary.h"
#include "timeloop.h"

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
/*
    __global__
    void coriolis_2nd_g(double* const __restrict__ ut, double* const __restrict__ vt,
                        double* const __restrict__ u,  double* const __restrict__ v,
                        double* const __restrict__ ug, double* const __restrict__ vg,
                        const double fc, const double ugrid, const double vgrid,
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

    __global__
    void coriolis_4th_g(double* const __restrict__ ut, double* const __restrict__ vt,
                        double* const __restrict__ u,  double* const __restrict__ v,
                        double* const __restrict__ ug, double* const __restrict__ vg,
                        const double fc, const double ugrid, const double vgrid,
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

    __global__
    void advec_wls_2nd_g(double* const __restrict__ st, double* const __restrict__ s,
                         const double* const __restrict__ wls, const double* const __restrict__ dzhi,
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

    __global__
    void large_scale_source_g(double* const __restrict__ st, double* const __restrict__ sls,
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

    __global__
    void nudging_tendency_g(double* const __restrict__ st, double* const __restrict__ smn,
			    double* const __restrict__ snudge, double* const __restrict__ nudge_fac,
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

    __global__
    void update_time_dependent_prof_g(double* const __restrict__ prof, const double* const __restrict__ data,
                                      const double fac0, const double fac1,
                                      const int index0,  const int index1,
                                      const int kmax,    const int kgc)
    {
        const int k = blockIdx.x*blockDim.x + threadIdx.x;
        const int kk = kmax;

        if (k < kmax)
            prof[k+kgc] = fac0*data[index0*kk+k] + fac1*data[index1*kk+k];
    }
*/
} // end namespace

template<typename TF>
void Force<TF>::prepare_device()
{
/*
    const int nmemsize  = gd.kcells*sizeof(TF);

    if (swlspres == "geo")
    {
        cuda_safe_call(cudaMalloc(&ug_g, nmemsize));
        cuda_safe_call(cudaMalloc(&vg_g, nmemsize));

        cuda_safe_call(cudaMemcpy(ug_g, ug, nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(vg_g, vg, nmemsize, cudaMemcpyHostToDevice));
        if (swtimedep_geo == "1")
        {

            for (std::map<std::string, double *>::const_iterator it=timedepdata_geo.begin(); it!=timedepdata_geo.end(); ++it)
            {
                int nmemsize2 = gd.kmax*timedeptime_geo[it->first].size()*sizeof(double);
                cuda_safe_call(cudaMalloc(&timedepdata_geo_g[it->first], nmemsize2));
                cuda_safe_call(cudaMemcpy(timedepdata_geo_g[it->first], timedepdata_geo[it->first], nmemsize2, cudaMemcpyHostToDevice));
            }
        }
    }

    if (swls == "1")
    {
        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
        {
            cuda_safe_call(cudaMalloc(&lsprofs_g[*it], nmemsize));
            cuda_safe_call(cudaMemcpy(lsprofs_g[*it], lsprofs[*it], nmemsize, cudaMemcpyHostToDevice));
        }
        if (swtimedep_ls == "1")
        {
            for (std::map<std::string, double *>::const_iterator it=timedepdata_ls.begin(); it!=timedepdata_ls.end(); ++it)
            {
                int nmemsize2 = gd.kmax*timedeptime_ls[it->first].size()*sizeof(double);
                cuda_safe_call(cudaMalloc(&timedepdata_ls_g[it->first], nmemsize2));
                cuda_safe_call(cudaMemcpy(timedepdata_ls_g[it->first], timedepdata_ls[it->first], nmemsize2, cudaMemcpyHostToDevice));
            }
        }
    }

    if (swnudge == "1")
    {
        for (std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
        {
            cuda_safe_call(cudaMalloc(&nudgeprofs_g[*it], nmemsize));
            cuda_safe_call(cudaMemcpy(nudgeprofs_g[*it], nudgeprofs[*it], nmemsize, cudaMemcpyHostToDevice));
        }
        cuda_safe_call(cudaMalloc(&nudge_factor_g, nmemsize));
        cuda_safe_call(cudaMemcpy(nudge_factor_g, nudge_factor, nmemsize, cudaMemcpyHostToDevice));
        if (swtimedep_nudge == "1")
        {
            for (std::map<std::string, double *>::const_iterator it=timedepdata_nudge.begin(); it!=timedepdata_nudge.end(); ++it)
            {
                int nmemsize2 = gd.kmax*timedeptime_nudge[it->first].size()*sizeof(double);
                cuda_safe_call(cudaMalloc(&timedepdata_nudge_g[it->first], nmemsize2));
                cuda_safe_call(cudaMemcpy(timedepdata_nudge_g[it->first], timedepdata_nudge[it->first], nmemsize2, cudaMemcpyHostToDevice));
            }
        }
    }

    if (swwls == "1")
    {
        cuda_safe_call(cudaMalloc(&wls_g, nmemsize));
        cuda_safe_call(cudaMemcpy(wls_g, wls, nmemsize, cudaMemcpyHostToDevice));
        if (swtimedep_wls == "1")
        {
            int nmemsize2 = gd.kmax*timedeptime_wls.size()*sizeof(double);
            cuda_safe_call(cudaMalloc(&timedepdata_wls_g, nmemsize2));
            cuda_safe_call(cudaMemcpy(timedepdata_wls_g, timedepdata_wls, nmemsize2, cudaMemcpyHostToDevice));
        }
    }
*/
}

template<typename TF>
void Force<TF>::clear_device()
{
/*    if (swlspres == "geo")
    {
        cuda_safe_call(cudaFree(ug_g));
        cuda_safe_call(cudaFree(vg_g));
        if (swtimedep_geo == "1")
        {
            for (std::map<std::string, double *>::const_iterator it=timedepdata_geo.begin(); it!=timedepdata_geo.end(); ++it)
                cuda_safe_call(cudaFree(timedepdata_geo_g[it->first]));
        }
    }

    if (swls == "1")
    {
        for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
            cuda_safe_call(cudaFree(lsprofs_g[*it]));
        if (swtimedep_ls == "1")
        {
            for (std::map<std::string, double *>::const_iterator it=timedepdata_ls.begin(); it!=timedepdata_ls.end(); ++it)
                cuda_safe_call(cudaFree(timedepdata_ls_g[it->first]));
        }
    }

    if (swnudge == "1")
    {
        for(std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
            cuda_safe_call(cudaFree(nudgeprofs_g[*it]));
        cuda_safe_call(cudaFree(nudge_factor_g));
        if (swtimedep_nudge == "1")
        {
            for (std::map<std::string, double *>::const_iterator it=timedepdata_nudge.begin(); it!=timedepdata_nudge.end(); ++it)
                cuda_safe_call(cudaFree(timedepdata_nudge_g[it->first]));
        }
    }

    if (swwls == "1")
    {
        cuda_safe_call(cudaFree(wls_g));
        if (swtimedep_wls == "1")
        {
            cuda_safe_call(cudaFree(timedepdata_wls_g));
        }
   }
*/
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
/*
        flux_step_1_g<TF><<<gridGPU, blockGPU>>>(
            &tmp->fld_g[offs], &fields.mp.at("u")->fld_g[offs],
            gd.dz_g,
            gd.icells, gd.ijcells,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend,    gd.jend,   gd.kend);
        cuda_check_error();
*/
        TF uavg  = field3d_operators.calc_mean(fields.mp.at("u")->fld_g);
        TF utavg = field3d_operators.calc_mean(fields.mt.at("u")->fld_g);
/*        flux_step_1_g<TF><<<gridGPU, blockGPU>>>(
            &tmp->fld_g[offs], &fields.mt.at("u")->fld_g[offs],
            gd.dz_g,
            gd.icells, gd.ijcells,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend,    gd.jend,   gd.kend);
        cuda_check_error();
*/
//        TF utavg  = field3d_operators.calc_mean(&tmp->fld_g[offs]);
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
/*    else if (swlspres == "geo")
    {
        if (gd.swspatialorder == "2")
        {
            coriolis_2nd_g<<<gridGPU, blockGPU>>>(
                &fields->ut->data_g[offs], &fields->vt->data_g[offs],
                &fields->u->data_g[offs],  &fields->v->data_g[offs],
                ug_g, vg_g, fc, grid.utrans, grid.vtrans,
                gd.icells, gd.ijcells,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend);
            cuda_check_error();
        }
        else if (gd.swspatialorder == "4")
        {
            coriolis_4th_g<<<gridGPU, blockGPU>>>(
                &fields->ut->data_g[offs], &fields->vt->data_g[offs],
                &fields->u->data_g[offs],  &fields->v->data_g[offs],
                ug_g, vg_g, fc, grid.utrans, grid.vtrans,
                gd.icells, gd.ijcells,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend);
            cuda_check_error();
        }
    }
*/
/*
    if (swls == "1")
    {
        for (std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
        {
            large_scale_source_g<<<gridGPU, blockGPU>>>(
                &fields->st[*it]->data_g[offs], lsprofs_g[*it],
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }

    if (swnudge == "1")
    {
        for(std::vector<std::string>::const_iterator it=nudgelist.begin(); it!=nudgelist.end(); ++it)
        {
            nudging_tendency_g<<<gridGPU, blockGPU>>>(
                &fields->at[*it]->data_g[offs],  fields->ap[*it]->datamean_g,
                nudgeprofs_g[*it], nudge_factor_g,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }

    if (swwls == "1")
    {
        for (FieldMap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
        {
            advec_wls_2nd_g<<<gridGPU, blockGPU>>>(
                &it->second->data_g[offs], fields->sp[it->first]->datamean_g, wls_g, gd.dzhi_g,
                gd.istart,  gd.jstart, gd.kstart,
                gd.iend,    gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }
*/
}
#endif
/*
#ifdef USECUDA
void Force::update_time_dependent_profs(std::map<std::string, double*>& profiles, std::map<std::string, double*> time_profiles,
                                        std::map<std::string, std::vector<double>> times, std::string suffix)
{
    const int blockk = 128;
    const int gridk  = gd.kmax/blockk + (gd.kmax%blockk > 0);

    // Loop over all profiles which might be time dependent
    for (auto& it : profiles)
    {
        std::string name = it.first + suffix;

        // Check if they have time dependent data
        if (time_profiles.find(name) != time_profiles.end())
        {
            // Get/calculate the interpolation indexes/factors
            int index0, index1;
            double fac0, fac1;

            model->timeloop->get_interpolation_factors(index0, index1, fac0, fac1, times[name]);

            // Calculate the new vertical profile
            update_time_dependent_prof_g<<<gridk, blockk>>>(
                it.second, time_profiles[name], fac0, fac1, index0, index1, gd.kmax, gd.kgc);
            cuda_check_error();
        }
    }

}
#endif
*/
/*
#ifdef USECUDA
void Force::update_time_dependent_prof(double* const prof, const double* const data, std::vector<double> times)
{
    const int blockk = 128;
    const int gridk  = gd.kmax/blockk + (gd.kmax%blockk > 0);

    int index0, index1;
    double fac0, fac1;

    model->timeloop->get_interpolation_factors(index0, index1, fac0, fac1, times);

    // Calculate the new vertical profile
    update_time_dependent_prof_g<<<gridk, blockk>>>(
        prof, data, fac0, fac1, index0, index1, gd.kmax, gd.kgc);
    cuda_check_error();

}
#endif
*/

template class Force<double>;
template class Force<float>;
