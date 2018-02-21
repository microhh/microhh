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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include <cufft.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres_2.h"
#include "defines.h"
#include "tools.h"

#ifdef USECUDA
template<typename TF>
void Pres_2<TF>::exec(const double dt)
{
/*    const Grid_data<TF>& gd = grid.get_grid_data();

    // create the input for the pressure solver
    input(fields.sd.at("p")->fld.data(),
          fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
          fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
          gd.dzi.data(), fields.rhoref.data(), fields.rhorefh.data(),
          dt);

    // solve the system
    solve(fields.sd.at("p")->fld.data(), fields.atmp.at("tmp1")->fld.data(), fields.atmp.at("tmp2")->fld.data(),
          gd.dz.data(), fields.rhoref.data(),
          grid.fftini, grid.fftouti, grid.fftinj, grid.fftoutj);

    // get the pressure tendencies from the pressure field
    output(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
           fields.sd.at("p")->fld.data(), gd.dzhi.data());
*/
}
#endif

#ifdef USECUDA
template<typename TF>
TF Pres_2<TF>::check_divergence()
{
/*    const Grid_data<TF>& gd = grid.get_grid_data();
    return calc_divergence(fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                           gd.dzi.data(), fields.rhoref.data(), fields.rhorefh.data());
*/
    return 0;
}
#endif


#ifdef USECUDA
template<typename TF>
void Pres_2<TF>::prepare_device()
{
/*    const int kmemsize = grid->kmax*sizeof(double);
    const int imemsize = grid->itot*sizeof(double);
    const int jmemsize = grid->jtot*sizeof(double);

    const int ijmemsize = grid->imax*grid->jmax*sizeof(double);

    cuda_safe_call(cudaMalloc((void**)&bmati_g, imemsize  ));
    cuda_safe_call(cudaMalloc((void**)&bmatj_g, jmemsize  ));
    cuda_safe_call(cudaMalloc((void**)&a_g, kmemsize      ));
    cuda_safe_call(cudaMalloc((void**)&c_g, kmemsize      ));
    cuda_safe_call(cudaMalloc((void**)&work2d_g, ijmemsize));

    cuda_safe_call(cudaMemcpy(bmati_g, bmati, imemsize, cudaMemcpyHostToDevice   ));
    cuda_safe_call(cudaMemcpy(bmatj_g, bmatj, jmemsize, cudaMemcpyHostToDevice   ));
    cuda_safe_call(cudaMemcpy(a_g, a, kmemsize, cudaMemcpyHostToDevice           ));
    cuda_safe_call(cudaMemcpy(c_g, c, kmemsize, cudaMemcpyHostToDevice           ));
    cuda_safe_call(cudaMemcpy(work2d_g, work2d, ijmemsize, cudaMemcpyHostToDevice));

    make_cufft_plan();
*/}
#endif

#ifdef USECUDA
template<typename TF>
void Pres_2<TF>::clear_device()
{
/*    cuda_safe_call(cudaFree(bmati_g ));
    cuda_safe_call(cudaFree(bmatj_g ));
    cuda_safe_call(cudaFree(a_g     ));
    cuda_safe_call(cudaFree(c_g     ));
    cuda_safe_call(cudaFree(work2d_g));
*/}
#endif

#ifdef USECUDA
template<typename TF>
TF Pres_2<TF>::calc_divergence(const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
                               const TF* const restrict dzi,
                               const TF* const restrict rhoref, const TF* const restrict rhorefh)
{
/*    const Grid_data<TF>& gd = grid.get_grid_data();

    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const TF dxi = 1./gd.dx;
    const TF dyi = 1./gd.dy;

    TF div = 0.;
    TF divmax = 0.;

    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                div = rhoref[k]*((u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi)
                    + (rhorefh[k+1]*w[ijk+kk]-rhorefh[k]*w[ijk])*dzi[k];

                divmax = std::max(divmax, std::abs(div));
            }

    master.max(&divmax, 1);

    return divmax;
*/  return 0;
}
#endif

template class Pres_2<double>;
template class Pres_2<float>;
