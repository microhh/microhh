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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "fft.h"
#include "pres_2.h"
#include "defines.h"
#include "stats.h"

template<typename TF>
Pres_2<TF>::Pres_2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, FFT<TF>& fftin, Input& inputin) :
    Pres<TF>(masterin, gridin, fieldsin, fftin, inputin),
    boundary_cyclic(master, grid)
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
    //clear_device();
    #endif
}

template<typename TF>
void Pres_2<TF>::create(Stats<TF>& stats)
{
    stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
void Pres_2<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // create the input for the pressure solver
    input(fields.sd.at("p")->fld.data(),
          fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
          fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
          gd.dzi.data(), fields.rhoref.data(), fields.rhorefh.data(),
          dt);

    // solve the system
    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    solve(fields.sd.at("p")->fld.data(), tmp1->fld.data(), tmp2->fld.data(),
          gd.dz.data(), fields.rhoref.data());

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);

    // get the pressure tendencies from the pressure field
    output(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
           fields.sd.at("p")->fld.data(), gd.dzhi.data());

   stats.calc_tend(*fields.mt.at("u"), tend_name);
   stats.calc_tend(*fields.mt.at("v"), tend_name);
   stats.calc_tend(*fields.mt.at("w"), tend_name);
}
#endif

#ifndef USECUDA
template<typename TF>
TF Pres_2<TF>::check_divergence()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    return calc_divergence(fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
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

    boundary_cyclic.init();
    fft.init();
}

template<typename TF>
void Pres_2<TF>::set_values()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    // Compute the modified wave numbers of the 2nd order scheme.
    const TF dxidxi = 1./(gd.dx*gd.dx);
    const TF dyidyi = 1./(gd.dy*gd.dy);

    const TF pi = std::acos(-1.);

    for (int j=0; j<gd.jtot/2+1; ++j)
        bmatj[j] = 2. * (std::cos(2.*pi*(TF)j/(TF)gd.jtot)-1.) * dyidyi;

    for (int j=gd.jtot/2+1; j<gd.jtot; ++j)
        bmatj[j] = bmatj[gd.jtot-j];

    for (int i=0; i<gd.itot/2+1; ++i)
        bmati[i] = 2. * (std::cos(2.*pi*(TF)i/(TF)gd.itot)-1.) * dxidxi;

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
void Pres_2<TF>::input(TF* const restrict p,
                       const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
                       TF* const restrict ut, TF* const restrict vt, TF* const restrict wt,
                       const TF* const restrict dzi, const TF* const restrict rhoref, const TF* const restrict rhorefh,
                       const TF dt)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const int jjp = gd.imax;
    const int kkp = gd.imax*gd.jmax;

    const TF dxi = TF(1.)/gd.dx;
    const TF dyi = TF(1.)/gd.dy;
    const TF dti = TF(1.)/dt;

    const int igc = gd.igc;
    const int jgc = gd.jgc;
    const int kgc = gd.kgc;

    // set the cyclic boundary conditions for the tendencies
    boundary_cyclic.exec(ut, Edge::East_west_edge  );
    boundary_cyclic.exec(vt, Edge::North_south_edge);

    // write pressure as a 3d array without ghost cells
    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                const int ijkp = i + j*jjp + k*kkp;
                const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
                p[ijkp] = rhoref[k+kgc] * ( (ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti) ) * dxi
                        + rhoref[k+kgc] * ( (vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti) ) * dyi
                        + ( rhorefh[k+kgc+1] * (wt[ijk+kk] + w[ijk+kk] * dti)
                          - rhorefh[k+kgc  ] * (wt[ijk   ] + w[ijk   ] * dti) ) * dzi[k+kgc];
            }
}

namespace
{
    // tridiagonal matrix solver, taken from Numerical Recipes, Press
    template<typename TF>
    void tdma(TF* const restrict a, TF* const restrict b, TF* const restrict c,
              TF* const restrict p, TF* const restrict work2d, TF* const restrict work3d,
              const int iblock, const int jblock, const int kmax)

    {
        const int jj = iblock;
        const int kk = iblock*jblock;

        for (int j=0; j<jblock; j++)
            #pragma ivdep
            for (int i=0; i<iblock; i++)
            {
                const int ij = i + j*jj;
                work2d[ij] = b[ij];
            }

        for (int j=0; j<jblock; j++)
            #pragma ivdep
            for (int i=0; i<iblock; i++)
            {
                const int ij = i + j*jj;
                p[ij] /= work2d[ij];
            }

        for (int k=1; k<kmax; k++)
        {
            for (int j=0; j<jblock; j++)
                #pragma ivdep
                for (int i=0; i<iblock; i++)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;
                    work3d[ijk] = c[k-1] / work2d[ij];
                }
            for (int j=0; j<jblock; j++)
                #pragma ivdep
                for (int i=0; i<iblock; i++)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;
                    work2d[ij] = b[ijk] - a[k]*work3d[ijk];
                }
            for (int j=0; j<jblock; j++)
                #pragma ivdep
                for (int i=0; i<iblock; i++)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;
                    p[ijk] -= a[k]*p[ijk-kk];
                    p[ijk] /= work2d[ij];
                }
        }

        for (int k=kmax-2; k>=0; k--)
            for (int j=0; j<jblock; j++)
                #pragma ivdep
                for (int i=0; i<iblock; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    p[ijk] -= work3d[ijk+kk]*p[ijk+kk];
                }
    }
}

template<typename TF>
void Pres_2<TF>::solve(TF* const restrict p, TF* const restrict work3d, TF* const restrict b,
                       const TF* const restrict dz, const TF* const restrict rhoref)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int imax   = gd.imax;
    const int jmax   = gd.jmax;
    const int kmax   = gd.kmax;
    const int iblock = gd.iblock;
    const int jblock = gd.jblock;
    const int igc    = gd.igc;
    const int jgc    = gd.jgc;
    const int kgc    = gd.kgc;

    int i,j,k,jj,kk,ijk;
    int iindex,jindex;

    fft.exec_forward(p, work3d);

    jj = iblock;
    kk = iblock*jblock;

    // solve the tridiagonal system
    // create vectors that go into the tridiagonal matrix solver
    for (k=0; k<kmax; k++)
        for (j=0; j<jblock; j++)
            #pragma ivdep
            for (i=0; i<iblock; i++)
            {
                // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
                iindex = md.mpicoordy * iblock + i;
                jindex = md.mpicoordx * jblock + j;

                ijk  = i + j*jj + k*kk;
                b[ijk] = dz[k+kgc]*dz[k+kgc] * rhoref[k+kgc]*(bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
                p[ijk] = dz[k+kgc]*dz[k+kgc] * p[ijk];
            }

    for (j=0; j<jblock; j++)
        #pragma ivdep
        for (i=0; i<iblock; i++)
        {
            iindex = md.mpicoordy * iblock + i;
            jindex = md.mpicoordx * jblock + j;

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
    tdma(a.data(), b, c.data(), p, work2d.data(), work3d,
         gd.iblock, gd.jblock, gd.kmax);

    fft.exec_backward(p, work3d);

    jj = imax;
    kk = imax*jmax;

    int ijkp,jjp,kkp;
    jjp = gd.icells;
    kkp = gd.ijcells;

    // put the pressure back onto the original grid including ghost cells
    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
                ijk  = i + j*jj + k*kk;
                p[ijkp] = work3d[ijk];
            }

    // set the boundary conditions
    // set a zero gradient boundary at the bottom
    for (int j=gd.jstart; j<gd.jend; ++j)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; ++i)
        {
            ijk = i + j*jjp + gd.kstart*kkp;
            p[ijk-kkp] = p[ijk];
        }

    // set the cyclic boundary conditions
    boundary_cyclic.exec(p);
}

template<typename TF>
void Pres_2<TF>::output(TF* const restrict ut, TF* const restrict vt, TF* const restrict wt,
                        const TF* const restrict p, const TF* const restrict dzhi)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const TF dxi = TF(1.)/gd.dx;
    const TF dyi = TF(1.)/gd.dy;

    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
                vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
                wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
            }
}

#ifndef USECUDA
template<typename TF>
TF Pres_2<TF>::calc_divergence(const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
                               const TF* const restrict dzi,
                               const TF* const restrict rhoref, const TF* const restrict rhorefh)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const TF dxi = TF(1.)/gd.dx;
    const TF dyi = TF(1.)/gd.dy;

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
}
#endif


#ifdef FLOAT_SINGLE
template class Pres_2<float>;
#else
template class Pres_2<double>;
#endif
