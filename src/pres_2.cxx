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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres_2.h"
#include "defines.h"
#include "model.h"

template<typename TF>
Pres_2<TF>::Pres_2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) : 
    Pres<TF>(masterin, gridin, fieldsin, inputin)
{
    #ifdef USECUDA
    a_g = 0;
    c_g = 0;
    work2d_g = 0;
    bmati_g  = 0;
    bmatj_g  = 0;
    #endif
}

template<typename TF>
Pres_2<TF>::~Pres_2()
{
    #ifdef USECUDA
    clear_device();
    #endif
}

#ifndef USECUDA
template<typename TF>
void Pres_2<TF>::exec(const double dt)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    // create the input for the pressure solver
    input(fields.sd.at("p")->data.data(),
          fields.mp.at("u")->data.data(), fields.mp.at("v")->data.data(), fields.mp.at("w")->data.data(),
          fields.mt.at("u")->data.data(), fields.mt.at("v")->data.data(), fields.mt.at("w")->data.data(),
          gd.dzi.data(), fields.rhoref.data(), fields.rhorefh.data(),
          dt);

    // solve the system
    solve(fields.sd.at("p")->data.data(), fields.atmp.at("tmp1")->data.data(), fields.atmp.at("tmp2")->data.data(),
          gd.dz.data(), fields.rhoref.data(),
          grid.fftini, grid.fftouti, grid.fftinj, grid.fftoutj);

    // get the pressure tendencies from the pressure field
    output(fields.mt.at("u")->data.data(), fields.mt.at("v")->data.data(), fields.mt.at("w")->data.data(), 
           fields.sd.at("p")->data.data(), gd.dzhi.data());
}
#endif

#ifndef USECUDA
template<typename TF>
double Pres_2<TF>::check_divergence()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    return calc_divergence(fields.mp.at("u")->data.data(), fields.mp.at("v")->data.data(), fields.mp.at("w")->data.data(),
                           gd.dzi.data(), fields.rhoref.data(), fields.rhorefh.data());
}
#endif

template<typename TF>
void Pres_2<TF>::init()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    bmati.resize(gd.itot);
    bmatj.resize(gd.jtot);

    a.resize(gd.kmax);
    c.resize(gd.kmax);

    work2d.resize(gd.imax*gd.jmax);
}

template<typename TF>
void Pres_2<TF>::set_values()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    // Compute the modified wave numbers of the 2nd order scheme.
    const double dxidxi = 1./(gd.dx*gd.dx);
    const double dyidyi = 1./(gd.dy*gd.dy);

    const double pi = std::acos(-1.);

    for (int j=0; j<gd.jtot/2+1; ++j)
        bmatj[j] = 2. * (std::cos(2.*pi*(double)j/(double)gd.jtot)-1.) * dyidyi;

    for (int j=gd.jtot/2+1; j<gd.jtot; ++j)
        bmatj[j] = bmatj[gd.jtot-j];

    for (int i=0; i<gd.itot/2+1; ++i)
        bmati[i] = 2. * (std::cos(2.*pi*(double)i/(double)gd.itot)-1.) * dxidxi;

    for (int i=gd.itot/2+1; i<gd.itot; ++i)
        bmati[i] = bmati[gd.itot-i];

    // create vectors that go into the tridiagonal matrix solver
    for (int k=0; k<gd.kmax; ++k)
    {
        a[k] = gd.dz[k+gd.kgc] * fields.rhorefh[k+gd.kgc  ]*gd.dzhi[k+gd.kgc  ];
        c[k] = gd.dz[k+gd.kgc] * fields.rhorefh[k+gd.kgc+1]*gd.dzhi[k+gd.kgc+1];
    }
}

template<typename TF>
void Pres_2<TF>::input(double* const restrict p, 
                       const double* const restrict u , const double* const restrict v , const double* const restrict w ,
                       const double* const restrict ut, const double* const restrict vt, const double* const restrict wt,
                       const double* const restrict dzi, const double* const restrict rhoref, const double* const restrict rhorefh,
                       const double dt)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const int jjp = grid->imax;
    const int kkp = grid->imax*grid->jmax;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;
    const double dti = 1./dt;

    const int igc = grid->igc;
    const int jgc = grid->jgc;
    const int kgc = grid->kgc;

    // set the cyclic boundary conditions for the tendencies
    grid->boundary_cyclic(ut, Edge::East_west_edge  );
    grid->boundary_cyclic(vt, Edge::North_south_edge);

    // write pressure as a 3d array without ghost cells
    for (int k=0; k<grid->kmax; k++)
        for (int j=0; j<grid->jmax; j++)
#pragma ivdep
            for (int i=0; i<grid->imax; i++)
            {
                const int ijkp = i + j*jjp + k*kkp;
                const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
                p[ijkp] = rhoref[k+kgc] * ( (ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti) ) * dxi
                        + rhoref[k+kgc] * ( (vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti) ) * dyi
                        + ( rhorefh[k+kgc+1] * (wt[ijk+kk] + w[ijk+kk] * dti) 
                          - rhorefh[k+kgc  ] * (wt[ijk   ] + w[ijk   ] * dti) ) * dzi[k+kgc];
            }
}

template<typename TF>
void Pres_2<TF>::solve(double* const restrict p, double* const restrict work3d, double* const restrict b,
                       const double* const restrict dz, const double* const restrict rhoref,
                       double* const restrict fftini, double* const restrict fftouti, 
                       double* const restrict fftinj, double* const restrict fftoutj)
{
    const int imax   = grid->imax;
    const int jmax   = grid->jmax;
    const int kmax   = grid->kmax;
    const int iblock = grid->iblock;
    const int jblock = grid->jblock;
    const int igc    = grid->igc;
    const int jgc    = grid->jgc;
    const int kgc    = grid->kgc;

    int i,j,k,jj,kk,ijk;
    int iindex,jindex;

    grid->fft_forward(p, work3d, fftini, fftouti, fftinj, fftoutj);

    jj = iblock;
    kk = iblock*jblock;

    //for (int i=0; i<itot*jtot; ++i)
    //  printf("%i %e\n",i,p[i]);
    //exit(1);

    // solve the tridiagonal system
    // create vectors that go into the tridiagonal matrix solver
    for (k=0; k<kmax; k++)
        for (j=0; j<jblock; j++)
#pragma ivdep
            for (i=0; i<iblock; i++)
            {
                // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
                iindex = master->mpicoordy * iblock + i;
                jindex = master->mpicoordx * jblock + j;

                ijk  = i + j*jj + k*kk;
                b[ijk] = dz[k+kgc]*dz[k+kgc] * rhoref[k+kgc]*(bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
                p[ijk] = dz[k+kgc]*dz[k+kgc] * p[ijk];
            }

    for (j=0; j<jblock; j++)
#pragma ivdep
        for (i=0; i<iblock; i++)
        {
            iindex = master->mpicoordy * iblock + i;
            jindex = master->mpicoordx * jblock + j;

            // substitute BC's
            ijk = i + j*jj;
            b[ijk] += a[0];

            // for wave number 0, which contains average, set pressure at top to zero
            ijk  = i + j*jj + (kmax-1)*kk;
            if (iindex == 0 && jindex == 0)
                b[ijk] -= c[kmax-1];
            // set dp/dz at top to zero
            else
                b[ijk] += c[kmax-1];
        }

    // call tdma solver
    tdma(a, b, c, p, work2d, work3d);

    grid->fft_backward(p, work3d, fftini, fftouti, fftinj, fftoutj);

    jj = imax;
    kk = imax*jmax;

    int ijkp,jjp,kkp;
    jjp = grid->icells;
    kkp = grid->ijcells;

    // put the pressure back onto the original grid including ghost cells
    for (int k=0; k<grid->kmax; k++)
        for (int j=0; j<grid->jmax; j++)
#pragma ivdep
            for (int i=0; i<grid->imax; i++)
            {
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
                ijk  = i + j*jj + k*kk;
                p[ijkp] = work3d[ijk];
            }

    // set the boundary conditions
    // set a zero gradient boundary at the bottom
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            ijk = i + j*jjp + grid->kstart*kkp;
            p[ijk-kkp] = p[ijk];
        }

    // set the cyclic boundary conditions
    grid->boundary_cyclic(p);
}

template<typename TF>
void Pres_2<TF>::output(double* const restrict ut, double* const restrict vt, double* const restrict wt, 
                        const double* const restrict p , const double* const restrict dzhi)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    for (int k=grid->kstart; k<grid->kend; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
                vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
                wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
            }
}

// tridiagonal matrix solver, taken from Numerical Recipes, Press
template<typename TF>
void Pres_2<TF>::tdma(double* restrict a, double* restrict b, double* restrict c, 
                      double* restrict p, double* restrict work2d, double* restrict work3d)

{
    int i,j,k,jj,kk,ijk,ij;
    int iblock,jblock,kmax;

    iblock = grid->iblock;
    jblock = grid->jblock;
    kmax = grid->kmax;

    jj = iblock;
    kk = iblock*jblock;

    for (j=0;j<jblock;j++)
#pragma ivdep
        for (i=0;i<iblock;i++)
        {
            ij = i + j*jj;
            work2d[ij] = b[ij];
        }

    for (j=0;j<jblock;j++)
#pragma ivdep
        for (i=0;i<iblock;i++)
        {
            ij = i + j*jj;
            p[ij] /= work2d[ij];
        }

    for (k=1; k<kmax; k++)
    {
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + k*kk;
                work3d[ijk] = c[k-1] / work2d[ij];
            }
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + k*kk;
                work2d[ij] = b[ijk] - a[k]*work3d[ijk];
            }
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + k*kk;
                p[ijk] -= a[k]*p[ijk-kk];
                p[ijk] /= work2d[ij];
            }
    }

    for (k=kmax-2; k>=0; k--)
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ijk = i + j*jj + k*kk;
                p[ijk] -= work3d[ijk+kk]*p[ijk+kk];
            }
}

#ifndef USECUDA
template<typename TF>
double Pres_2<TF>::calc_divergence(const double* const restrict u, const double* const restrict v, const double* const restrict w,
                                   const double* const restrict dzi,
                                   const double* const restrict rhoref, const double* const restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double div    = 0.;
    double divmax = 0.;

    for (int k=grid->kstart; k<grid->kend; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                div = rhoref[k]*((u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi) 
                    + (rhorefh[k+1]*w[ijk+kk]-rhorefh[k]*w[ijk])*dzi[k];

                divmax = std::max(divmax, std::abs(div));
            }

    grid->get_max(&divmax);

    return divmax;
}
#endif

template class Pres_2<double>;
template class Pres_2<float>;
