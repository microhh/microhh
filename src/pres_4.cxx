/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
 *
 * The heptadiagonal matrix solver is
 * Copyright (c) 2014 Juan Pedro Mellado
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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres_4.h"
#include "defines.h"
#include "finite_difference.h"
#include "model.h"
#include "stats.h"

using namespace Finite_difference::O4;

template<typename TF>
Pres_4<TF>::Pres_4(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, FFT<TF>& fftin, Input& inputin) :
    Pres<TF>(masterin, gridin, fieldsin, fftin, inputin),
    boundary_cyclic(master, grid)
{
    #ifdef USECUDA
    bmati_g = 0;
    bmatj_g = 0;
    m1_g = 0;
    m2_g = 0;
    m3_g = 0;
    m4_g = 0;
    m5_g = 0;
    m6_g = 0;
    m7_g = 0;
    #endif
}

template<typename TF>
Pres_4<TF>::~Pres_4()
{
    #ifdef USECUDA
    // clear_device();
    #endif
}

template<typename TF>
void Pres_4<TF>::create(Stats<TF>& stats)
{
    stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
void Pres_4<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // 1. Create the input for the pressure solver.
    // In case of a two-dimensional run, remove calculation of v contribution.
    if (gd.jtot == 1)
        input<false>(fields.sd.at("p")->fld.data(),
                     fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                     fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
                     gd.dzi4.data(), dt);
    else
        input<true>(fields.sd.at("p")->fld.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                    fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
                    gd.dzi4.data(), dt);

    // 2. Solve the Poisson equation using FFTs and a heptadiagonal solver

    /* Find the thickness of a vectorizable slice. There is a need for 8 slices for the pressure
       solver and we use two three dimensional temp fields, so there are 4 slices per field.
       The thickness is therefore jblock/4. Since there are always three ghost cells, even in a 2D
       run the fields are large enough. */
    // const int jslice = std::max(grid->jblock/4, 1);

    /* The CPU version gives the best performance in case jslice = 1, due to cache misses.
       In case this value will be set to larger than 1, checks need to be build in for out of bounds
       reads in case jblock does not divide by 4. */
    const int jslice = 1;

    auto tmp1 = fields.get_tmp();

    auto tmp2_fld = fields.get_tmp();
    auto tmp3_fld = fields.get_tmp();

    // Shortcuts for simpler notation.
    TF* tmp2 = tmp2_fld->fld.data();
    TF* tmp3 = tmp3_fld->fld.data();

    const int ns = gd.iblock*jslice*(gd.kmax+4);

    solve(fields.sd.at("p")->fld.data(), tmp1->fld.data(), gd.dz.data(),
          m1.data(), m2.data(), m3.data(), m4.data(),
          m5.data(), m6.data(), m7.data(),
          &tmp2[0*ns], &tmp2[1*ns], &tmp2[2*ns], &tmp2[3*ns],
          &tmp3[0*ns], &tmp3[1*ns], &tmp3[2*ns], &tmp3[3*ns],
          bmati.data(), bmatj.data(),
          jslice);

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2_fld);
    fields.release_tmp(tmp3_fld);

    // 3. Get the pressure tendencies from the pressure field.
    if (gd.jtot == 1)
        output<false>(
                fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
                fields.sd.at("p")->fld.data(), gd.dzhi4.data());
    else
        output<true>(
                fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
                fields.sd.at("p")->fld.data(), gd.dzhi4.data());

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);

}

template<typename TF>
TF Pres_4<TF>::check_divergence()
{
    auto& gd = grid.get_grid_data();
    return calc_divergence(
                    fields.mp.at("u")->fld.data(),
                    fields.mp.at("v")->fld.data(),
                    fields.mp.at("w")->fld.data(),
                    gd.dzi4.data());
}
#endif

template<typename TF>
void Pres_4<TF>::init()
{
    auto& gd = grid.get_grid_data();

    bmati.resize(gd.itot);
    bmatj.resize(gd.jtot);

    m1.resize(gd.kmax);
    m2.resize(gd.kmax);
    m3.resize(gd.kmax);
    m4.resize(gd.kmax);
    m5.resize(gd.kmax);
    m6.resize(gd.kmax);
    m7.resize(gd.kmax);

    boundary_cyclic.init();
    fft.init();
}

template<typename TF>
void Pres_4<TF>::set_values()
{
    auto& gd = grid.get_grid_data();

    const int itot = gd.itot;
    const int jtot = gd.jtot;
    const int kmax = gd.kmax;
    const int kstart = gd.kstart;

    // compute the modified wave numbers of the 4th order scheme
    const TF dxidxi = 1./(gd.dx*gd.dx);
    const TF dyidyi = 1./(gd.dy*gd.dy);

    const TF pi = std::acos(-1.);

    // Convert the coefficients to float after calculation.
    for (int j=0; j<jtot/2+1; j++)
        bmatj[j] = ( 2.* (1./576.)    * std::cos(6.*pi*(double)j/(double)jtot)
                   - 2.* (54./576.)   * std::cos(4.*pi*(double)j/(double)jtot)
                   + 2.* (783./576.)  * std::cos(2.*pi*(double)j/(double)jtot)
                   -     (1460./576.) ) * dyidyi;

    for (int j=jtot/2+1; j<jtot; j++)
        bmatj[j] = bmatj[jtot-j];

    for (int i=0; i<itot/2+1; i++)
        bmati[i] = ( 2.* (1./576.)    * std::cos(6.*pi*(double)i/(double)itot)
                   - 2.* (54./576.)   * std::cos(4.*pi*(double)i/(double)itot)
                   + 2.* (783./576.)  * std::cos(2.*pi*(double)i/(double)itot)
                   -     (1460./576.) ) * dxidxi;

    for (int i=itot/2+1; i<itot; i++)
        bmati[i] = bmati[itot-i];

    // Shortcuts for easier notation.
    auto& dzi4  = gd.dzi4;
    auto& dzhi4 = gd.dzhi4;

    int k,kc;
    // create vectors that go into the matrix solver
    // bottom boundary, taking into account that w is mirrored over the wall to conserve global momentum
    k  = 0;
    kc = kstart+k;
    m1[k] = 0.;
    m2[k] = (1./576.) * (                 -  27.*dzhi4[kc]                                      ) * dzi4[kc];
    m3[k] = (1./576.) * ( -1.*dzhi4[kc+1] + 729.*dzhi4[kc] +  27.*dzhi4[kc+1]                   ) * dzi4[kc];
    m4[k] = (1./576.) * ( 27.*dzhi4[kc+1] - 729.*dzhi4[kc] - 729.*dzhi4[kc+1] -  1.*dzhi4[kc+2] ) * dzi4[kc];
    m5[k] = (1./576.) * (-27.*dzhi4[kc+1] +  27.*dzhi4[kc] + 729.*dzhi4[kc+1] + 27.*dzhi4[kc+2] ) * dzi4[kc];
    m6[k] = (1./576.) * (  1.*dzhi4[kc+1]                  -  27.*dzhi4[kc+1] - 27.*dzhi4[kc+2] ) * dzi4[kc];
    m7[k] = (1./576.) * (                                                     +  1.*dzhi4[kc+2] ) * dzi4[kc];

    for (int k=1; k<kmax-1; k++)
    {
        kc = kstart+k;
        m1[k] = (1./576.) * (   1.*dzhi4[kc-1]                                                       ) * dzi4[kc];
        m2[k] = (1./576.) * ( -27.*dzhi4[kc-1] -  27.*dzhi4[kc]                                      ) * dzi4[kc];
        m3[k] = (1./576.) * (  27.*dzhi4[kc-1] + 729.*dzhi4[kc] +  27.*dzhi4[kc+1]                   ) * dzi4[kc];
        m4[k] = (1./576.) * (  -1.*dzhi4[kc-1] - 729.*dzhi4[kc] - 729.*dzhi4[kc+1] -  1.*dzhi4[kc+2] ) * dzi4[kc];
        m5[k] = (1./576.) * (                  +  27.*dzhi4[kc] + 729.*dzhi4[kc+1] + 27.*dzhi4[kc+2] ) * dzi4[kc];
        m6[k] = (1./576.) * (                                   -  27.*dzhi4[kc+1] - 27.*dzhi4[kc+2] ) * dzi4[kc];
        m7[k] = (1./576.) * (                                                      +  1.*dzhi4[kc+2] ) * dzi4[kc];
    }

    // top boundary, taking into account that w is mirrored over the wall to conserve global momentum
    k  = kmax-1;
    kc = kstart+k;
    m1[k] = (1./576.) * (   1.*dzhi4[kc-1]                                                     ) * dzi4[kc];
    m2[k] = (1./576.) * ( -27.*dzhi4[kc-1] -  27.*dzhi4[kc]                    +  1.*dzhi4[kc] ) * dzi4[kc];
    m3[k] = (1./576.) * (  27.*dzhi4[kc-1] + 729.*dzhi4[kc] +  27.*dzhi4[kc+1] - 27.*dzhi4[kc] ) * dzi4[kc];
    m4[k] = (1./576.) * (  -1.*dzhi4[kc-1] - 729.*dzhi4[kc] - 729.*dzhi4[kc+1] + 27.*dzhi4[kc] ) * dzi4[kc];
    m5[k] = (1./576.) * (                  +  27.*dzhi4[kc] + 729.*dzhi4[kc+1] -  1.*dzhi4[kc] ) * dzi4[kc];
    m6[k] = (1./576.) * (                                   -  27.*dzhi4[kc+1]                 ) * dzi4[kc];
    m7[k] = 0.;
}

template<typename TF>
template<bool dim3>
void Pres_4<TF>::input(
        TF* restrict p,
        const TF* restrict u, const TF* restrict v, const TF* restrict w ,
        TF* restrict ut, TF* restrict vt, TF* restrict wt,
        const TF* restrict dzi4, const TF dt)
{
    auto& gd = grid.get_grid_data();

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*gd.icells;
    const int jj2 = 2*gd.icells;
    const int kk1 = 1*gd.ijcells;
    const int kk2 = 2*gd.ijcells;

    const int jjp = gd.imax;
    const int kkp = gd.imax*gd.jmax;

    const TF dxi = 1./gd.dx;
    const TF dyi = 1./gd.dy;
    const TF dti = 1./dt;

    const int igc = gd.igc;
    const int jgc = gd.jgc;
    const int kgc = gd.kgc;

    const int kmax = gd.kmax;

    // Set the cyclic boundary conditions for the tendencies.
    boundary_cyclic.exec(ut, Edge::East_west_edge);
    if (dim3)
        boundary_cyclic.exec(vt, Edge::North_south_edge);

    // Set the bc.
    for (int j=0; j<gd.jmax; j++)
        #pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            const int ijk  = i+igc + (j+jgc)*jj1 + kgc*kk1;
            wt[ijk-kk1] = -wt[ijk+kk1];
        }
    for (int j=0; j<gd.jmax; j++)
        #pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            const int ijk  = i+igc + (j+jgc)*jj1 + (kmax+kgc)*kk1;
            wt[ijk+kk1] = -wt[ijk-kk1];
        }

    for (int k=0; k<gd.kmax; k++)
        for (int j=0; j<gd.jmax; j++)
            #pragma ivdep
            for (int i=0; i<gd.imax; i++)
            {
                const int ijkp = i + j*jjp + k*kkp;
                const int ijk  = i+igc + (j+jgc)*jj1 + (k+kgc)*kk1;
                p[ijkp]  = (cg0<TF>*(ut[ijk-ii1] + u[ijk-ii1]*dti) + cg1<TF>*(ut[ijk] + u[ijk]*dti) + cg2<TF>*(ut[ijk+ii1] + u[ijk+ii1]*dti) + cg3<TF>*(ut[ijk+ii2] + u[ijk+ii2]*dti)) * dxi;
                if (dim3)
                    p[ijkp] += (cg0<TF>*(vt[ijk-jj1] + v[ijk-jj1]*dti) + cg1<TF>*(vt[ijk] + v[ijk]*dti) + cg2<TF>*(vt[ijk+jj1] + v[ijk+jj1]*dti) + cg3<TF>*(vt[ijk+jj2] + v[ijk+jj2]*dti)) * dyi;
                p[ijkp] += (cg0<TF>*(wt[ijk-kk1] + w[ijk-kk1]*dti) + cg1<TF>*(wt[ijk] + w[ijk]*dti) + cg2<TF>*(wt[ijk+kk1] + w[ijk+kk1]*dti) + cg3<TF>*(wt[ijk+kk2] + w[ijk+kk2]*dti)) * dzi4[k+kgc];
            }
}

template<typename TF>
void Pres_4<TF>::solve(
        TF* restrict p, TF* restrict work3d, const TF* restrict dz,
        const TF* restrict m1, const TF* restrict m2, const TF* restrict m3, const TF* restrict m4,
        const TF* restrict m5, const TF* restrict m6, const TF* restrict m7,
        TF* restrict m1temp, TF* restrict m2temp, TF* restrict m3temp, TF* restrict m4temp,
        TF* restrict m5temp, TF* restrict m6temp, TF* restrict m7temp, TF* restrict ptemp,
        TF* restrict bmati, TF* restrict bmatj,
        const int jslice)
{
    auto& gd = grid.get_grid_data();

    const int imax   = gd.imax;
    const int jmax   = gd.jmax;
    const int kmax   = gd.kmax;
    const int iblock = gd.iblock;
    const int jblock = gd.jblock;
    const int igc    = gd.igc;
    const int jgc    = gd.jgc;
    const int kgc    = gd.kgc;

    fft.exec_forward(p, work3d);

    int jj,kk,ik,ijk;
    int iindex,jindex;

    jj = iblock;
    kk = iblock*jblock;

    auto& md = master.get_MPI_data();

    const int mpicoordx = md.mpicoordx;
    const int mpicoordy = md.mpicoordy;

    // Calculate the step size.
    const int nj = jblock/jslice;

    const int kki1 = 1*iblock*jslice;
    const int kki2 = 2*iblock*jslice;
    const int kki3 = 3*iblock*jslice;

    for (int n=0; n<nj; ++n)
    {
        for (int j=0; j<jslice; ++j)
            #pragma ivdep
            for (int i=0; i<iblock; ++i)
            {
                // Set a zero gradient bc at the bottom.
                ik = i + j*jj;
                m1temp[ik] = TF( 0.);
                m2temp[ik] = TF( 0.);
                m3temp[ik] = TF( 0.);
                m4temp[ik] = TF( 1.);
                m5temp[ik] = TF( 0.);
                m6temp[ik] = TF( 0.);
                m7temp[ik] = TF(-1.);
                ptemp [ik] = TF( 0.);
            }

        for (int j=0; j<jslice; ++j)
            #pragma ivdep
            for (int i=0; i<iblock; ++i)
            {
                ik = i + j*jj;
                m1temp[ik+kki1] = TF( 0.);
                m2temp[ik+kki1] = TF( 0.);
                m3temp[ik+kki1] = TF( 0.);
                m4temp[ik+kki1] = TF( 1.);
                m5temp[ik+kki1] = TF(-1.);
                m6temp[ik+kki1] = TF( 0.);
                m7temp[ik+kki1] = TF( 0.);
                ptemp [ik+kki1] = TF( 0.);
            }

        for (int k=0; k<kmax; ++k)
            for (int j=0; j<jslice; ++j)
            {
                jindex = mpicoordx*jblock + n*jslice + j;
                #pragma ivdep
                for (int i=0; i<iblock; ++i)
                {
                    // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
                    iindex = mpicoordy*iblock + i;

                    ijk = i + (j + n*jslice)*jj + k*kk;
                    ik  = i + j*jj + k*kki1;
                    m1temp[ik+kki2] = m1[k];
                    m2temp[ik+kki2] = m2[k];
                    m3temp[ik+kki2] = m3[k];
                    m4temp[ik+kki2] = m4[k] + bmati[iindex] + bmatj[jindex];
                    m5temp[ik+kki2] = m5[k];
                    m6temp[ik+kki2] = m6[k];
                    m7temp[ik+kki2] = m7[k];
                    ptemp [ik+kki2] = p[ijk];
                }
            }

        for (int j=0; j<jslice; ++j)
        {
            jindex = mpicoordx*jblock + n*jslice + j;
            #pragma ivdep
            for (int i=0; i<iblock; ++i)
            {
                // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
                iindex = mpicoordy*iblock + i;

                // Set the top boundary.
                ik = i + j*jj + kmax*kki1;
                if (iindex == 0 && jindex == 0)
                {
                    m1temp[ik+kki2] = TF(   0.);
                    m2temp[ik+kki2] = TF(-1/3.);
                    m3temp[ik+kki2] = TF(   2.);
                    m4temp[ik+kki2] = TF(   1.);

                    m1temp[ik+kki3] = TF(  -2.);
                    m2temp[ik+kki3] = TF(   9.);
                    m3temp[ik+kki3] = TF(   0.);
                    m4temp[ik+kki3] = TF(   1.);
                }
                // Set dp/dz at top to zero.
                else
                {
                    m1temp[ik+kki2] = TF( 0.);
                    m2temp[ik+kki2] = TF( 0.);
                    m3temp[ik+kki2] = TF(-1.);
                    m4temp[ik+kki2] = TF( 1.);

                    m1temp[ik+kki3] = TF(-1.);
                    m2temp[ik+kki3] = TF( 0.);
                    m3temp[ik+kki3] = TF( 0.);
                    m4temp[ik+kki3] = TF( 1.);
                }
            }
        }

        for (int j=0; j<jslice; ++j)
            #pragma ivdep
            for (int i=0; i<iblock; ++i)
            {
                // Set the top boundary.
                ik = i + j*jj + kmax*kki1;
                m5temp[ik+kki2] = TF(0.);
                m6temp[ik+kki2] = TF(0.);
                m7temp[ik+kki2] = TF(0.);
                ptemp [ik+kki2] = TF(0.);

                m5temp[ik+kki3] = TF(0.);
                m6temp[ik+kki3] = TF(0.);
                m7temp[ik+kki3] = TF(0.);
                ptemp [ik+kki3] = TF(0.);
            }

        hdma(m1temp, m2temp, m3temp, m4temp, m5temp, m6temp, m7temp, ptemp, jslice);

        // Put back the solution.
        for (int k=0; k<kmax; ++k)
            for (int j=0; j<jslice; ++j)
                #pragma ivdep
                for (int i=0; i<iblock; ++i)
                {
                    const int ik  = i + j*jj + k*kki1;
                    const int ijk = i + (j + n*jslice)*jj + k*kk;
                    p[ijk] = ptemp[ik+kki2];
                }
    }

    fft.exec_backward(p, work3d);

    // Put the pressure back onto the original grid including ghost cells.
    jj = imax;
    kk = imax*jmax;

    int ijkp,jjp,kkp1,kkp2;
    jjp = gd.icells;
    kkp1 = 1*gd.ijcells;
    kkp2 = 2*gd.ijcells;

    for (int k=0; k<gd.kmax; k++)
        for (int j=0; j<gd.jmax; j++)
            #pragma ivdep
            for (int i=0; i<gd.imax; i++)
            {
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp1;
                ijk  = i + j*jj + k*kk;
                p[ijkp] = work3d[ijk];
            }

    // Set a zero gradient boundary at the bottom.
    for (int j=gd.jstart; j<gd.jend; j++)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; i++)
        {
            ijk = i + j*jjp + gd.kstart*kkp1;
            p[ijk-kkp1] = p[ijk     ];
            p[ijk-kkp2] = p[ijk+kkp1];
        }

    // Set a zero gradient boundary at the top.
    for (int j=gd.jstart; j<gd.jend; j++)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; i++)
        {
            ijk = i + j*jjp + (gd.kend-1)*kkp1;
            p[ijk+kkp1] = p[ijk     ];
            p[ijk+kkp2] = p[ijk-kkp1];
        }

    // Set the cyclic boundary conditions.
    boundary_cyclic.exec(p);
}

template<typename TF>
template<bool dim3>
void Pres_4<TF>::output(TF* restrict ut, TF* restrict vt, TF* restrict wt,
                        const TF* restrict p , const TF* restrict dzhi4)
{
    auto& gd = grid.get_grid_data();

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*gd.icells;
    const int jj2 = 2*gd.icells;
    const int kk1 = 1*gd.ijcells;
    const int kk2 = 2*gd.ijcells;

    const int kstart = gd.kstart;

    const TF dxi = 1./gd.dx;
    const TF dyi = 1./gd.dy;

    for (int j=gd.jstart; j<gd.jend; j++)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; i++)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            ut[ijk] -= (cg0<TF>*p[ijk-ii2] + cg1<TF>*p[ijk-ii1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+ii1]) * dxi;
            if (dim3)
                vt[ijk] -= (cg0<TF>*p[ijk-jj2] + cg1<TF>*p[ijk-jj1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+jj1]) * dyi;
        }

    for (int k=gd.kstart+1; k<gd.kend; k++)
        for (int j=gd.jstart; j<gd.jend; j++)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] -= (cg0<TF>*p[ijk-ii2] + cg1<TF>*p[ijk-ii1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+ii1]) * dxi;
                if (dim3)
                    vt[ijk] -= (cg0<TF>*p[ijk-jj2] + cg1<TF>*p[ijk-jj1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+jj1]) * dyi;
                wt[ijk] -= (cg0<TF>*p[ijk-kk2] + cg1<TF>*p[ijk-kk1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+kk1]) * dzhi4[k];
            }
}

template<typename TF>
void Pres_4<TF>::hdma(
        TF* restrict m1, TF* restrict m2, TF* restrict m3, TF* restrict m4,
        TF* restrict m5, TF* restrict m6, TF* restrict m7, TF* restrict p,
        const int jslice)
{
    auto& gd = grid.get_grid_data();

    const int kmax   = gd.kmax;
    const int iblock = gd.iblock;

    const int jj = gd.iblock;

    const int kk1 = 1*gd.iblock*jslice;
    const int kk2 = 2*gd.iblock*jslice;
    const int kk3 = 3*gd.iblock*jslice;

    int k,ik;

    // Use LU factorization.
    k = 0;
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj;
            m1[ik] = TF(1.);
            m2[ik] = TF(1.);
            m3[ik] = TF(1.)            / m4[ik];
            m4[ik] = TF(1.);
            m5[ik] = m5[ik]*m3[ik];
            m6[ik] = m6[ik]*m3[ik];
            m7[ik] = m7[ik]*m3[ik];
        }

    k = 1;
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj + k*kk1;
            m1[ik] = TF(1.);
            m2[ik] = TF(1.);
            m3[ik] = m3[ik]                     / m4[ik-kk1];
            m4[ik] = m4[ik] - m3[ik]*m5[ik-kk1];
            m5[ik] = m5[ik] - m3[ik]*m6[ik-kk1];
            m6[ik] = m6[ik] - m3[ik]*m7[ik-kk1];
        }

    k = 2;
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj + k*kk1;
            m1[ik] = TF(1.);
            m2[ik] =   m2[ik]                                           / m4[ik-kk2];
            m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] ) / m4[ik-kk1];
            m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2];
            m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
            m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
        }

    for (k=3; k<kmax+2; ++k)
        for (int j=0; j<jslice; ++j)
            #pragma ivdep
            for (int i=0; i<iblock; ++i)
            {
                ik = i + j*jj + k*kk1;
                m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
                m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
                m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
                m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
                m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
                m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
            }

    k = kmax+1;
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj + k*kk1;
            m7[ik] = TF(1.);
        }

    k = kmax+2;
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj + k*kk1;
            m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
            m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
            m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
            m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
            m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
            m6[ik] = TF(1.);
            m7[ik] = TF(1.);
        }

    k = kmax+3;
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj + k*kk1;
            m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
            m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
            m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
            m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
            m5[ik] = (1.);
            m6[ik] = (1.);
            m7[ik] = (1.);
        }

    // Do the backward substitution.
    // First, solve Ly = p, forward.
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj;
            p[ik    ] =             p[ik    ]*m3[ik    ];
            p[ik+kk1] = p[ik+kk1] - p[ik    ]*m3[ik+kk1];
            p[ik+kk2] = p[ik+kk2] - p[ik+kk1]*m3[ik+kk2] - p[ik]*m2[ik+kk2];
        }

    for (k=3; k<kmax+4; ++k)
        for (int j=0; j<jslice; ++j)
            #pragma ivdep
            for (int i=0; i<iblock; ++i)
            {
                ik = i + j*jj + k*kk1;
                p[ik] = p[ik] - p[ik-kk1]*m3[ik] - p[ik-kk2]*m2[ik] - p[ik-kk3]*m1[ik];
            }

    // Second, solve Ux=y, backward.
    k = kmax+3;
    for (int j=0; j<jslice; ++j)
        #pragma ivdep
        for (int i=0; i<iblock; ++i)
        {
            ik = i + j*jj + k*kk1;
            p[ik    ] =   p[ik    ]                                             / m4[ik    ];
            p[ik-kk1] = ( p[ik-kk1] - p[ik    ]*m5[ik-kk1] )                    / m4[ik-kk1];
            p[ik-kk2] = ( p[ik-kk2] - p[ik-kk1]*m5[ik-kk2] - p[ik]*m6[ik-kk2] ) / m4[ik-kk2];
        }

    for (k=kmax; k>=0; --k)
        for (int j=0; j<jslice; ++j)
            #pragma ivdep
            for (int i=0; i<iblock; ++i)
            {
                ik = i + j*jj + k*kk1;
                p[ik] = ( p[ik] - p[ik+kk1]*m5[ik] - p[ik+kk2]*m6[ik] - p[ik+kk3]*m7[ik] ) / m4[ik];
            }
}

template<typename TF>
TF Pres_4<TF>::calc_divergence(
        const TF* restrict u, const TF* restrict v, const TF* restrict w, const TF* restrict dzi4)
{
    auto& gd = grid.get_grid_data();

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*gd.icells;
    const int jj2 = 2*gd.icells;
    const int kk1 = 1*gd.ijcells;
    const int kk2 = 2*gd.ijcells;

    const TF dxi = 1./gd.dx;
    const TF dyi = 1./gd.dy;

    TF div, divmax;
    divmax = 0;

    for (int k=gd.kstart; k<gd.kend; k++)
        for (int j=gd.jstart; j<gd.jend; j++)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                div = (cg0<TF>*u[ijk-ii1] + cg1<TF>*u[ijk] + cg2<TF>*u[ijk+ii1] + cg3<TF>*u[ijk+ii2]) * dxi
                    + (cg0<TF>*v[ijk-jj1] + cg1<TF>*v[ijk] + cg2<TF>*v[ijk+jj1] + cg3<TF>*v[ijk+jj2]) * dyi
                    + (cg0<TF>*w[ijk-kk1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+kk1] + cg3<TF>*w[ijk+kk2]) * dzi4[k];

                divmax = std::max(divmax, std::abs(div));
            }

    master.max(&divmax, 1);

    return divmax;
}


#ifdef FLOAT_SINGLE
template class Pres_4<float>;
#else
template class Pres_4<double>;
#endif
