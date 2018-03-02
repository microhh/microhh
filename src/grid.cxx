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
#include <cstdlib>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "input.h"
#include "data_block.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"

/**
 * This function constructs the grid class.
 * @param modelin Pointer to the model class.
 * @param inputin Pointer to the input class.
 */
template<typename TF>
Grid<TF>::Grid(Master& masterin, Input& input) :
    master(masterin),
    transpose(master, *this)
{
    mpitypes = false;

    gd.xsize = input.get_item<TF>("grid", "xsize", "");
    gd.ysize = input.get_item<TF>("grid", "ysize", "");
    gd.zsize = input.get_item<TF>("grid", "zsize", "");

    gd.itot = input.get_item<int>("grid", "itot", "");
    gd.jtot = input.get_item<int>("grid", "jtot", "");
    gd.ktot = input.get_item<int>("grid", "ktot", "");

    utrans = input.get_item<TF>("grid", "utrans", "", 0.);
    vtrans = input.get_item<TF>("grid", "vtrans", "", 0.);

    swspatialorder = input.get_item<std::string>("grid", "swspatialorder", "");

    if (!(swspatialorder == "2" || swspatialorder == "4"))
    {
        master.print_error("\"%s\" is an illegal value for swspatialorder\n", swspatialorder.c_str());
        throw 1;
    }
    // 2nd order scheme requires only 1 ghost cell
    if (swspatialorder == "2")
    {
        gd.igc = 1;
        gd.jgc = 1;
        gd.kgc = 1;
    }
    // 4th order scheme requires 3 ghost cells
    else if (swspatialorder == "4")
    {
        gd.igc = 3;
        gd.jgc = 3;
        gd.kgc = 3;
    }
}

template<typename TF>
Grid<TF>::~Grid()
{
    exit_mpi();
}

/**
 * This function allocates the dynamic arrays in the field class
 * variables and calculates the derived grid indices and dimensions.
 */
template<typename TF>
void Grid<TF>::init()
{
    // Check whether the grid fits the processor configuration.
    if (gd.itot % master.npx != 0)
    {
        master.print_error("itot = %d is not a multiple of npx = %d\n", gd.itot, master.npx);
        throw 1;
    }
    if (gd.itot % master.npy != 0)
    {
        master.print_error("itot = %d is not a multiple of npy = %d\n", gd.itot, master.npy);
        throw 1;
    }
    // Check this one only when npy > 1, since the transpose in that direction only happens then.
    if (gd.jtot % master.npx != 0)
    {
        master.print_error("jtot = %d is not a multiple of npx = %d\n", gd.jtot, master.npx);
        throw 1;
    }
    if (gd.jtot % master.npy != 0)
    {
        master.print_error("jtot = %d is not a multiple of npy = %d\n", gd.jtot, master.npy);
        throw 1;
    }
    if (gd.ktot % master.npx != 0)
    {
        master.print_error("ERROR ktot = %d is not a multiple of npx = %d\n", gd.ktot, master.npx);
        throw 1;
    }

    // Calculate the total number of grid cells.
    gd.ntot = gd.itot*gd.jtot*gd.ktot;

    // Calculate the grid dimensions per process.
    gd.imax = gd.itot / master.npx;
    gd.jmax = gd.jtot / master.npy;
    gd.kmax = gd.ktot;
    gd.nmax = gd.imax*gd.jmax*gd.kmax;

    // Calculate the block sizes for the transposes.
    gd.iblock = gd.itot / master.npy;
    gd.jblock = gd.jtot / master.npx;
    gd.kblock = gd.ktot / master.npx;

    // Calculate the grid dimensions including ghost cells.
    gd.icells  = (gd.imax+2*gd.igc);
    gd.jcells  = (gd.jmax+2*gd.jgc);
    gd.ijcells = gd.icells*gd.jcells;
    gd.kcells  = (gd.kmax+2*gd.kgc);
    gd.ncells  = (gd.imax+2*gd.igc)*(gd.jmax+2*gd.jgc)*(gd.kmax+2*gd.kgc);

    // Calculate the starting and ending points for loops over the grid.
    gd.istart = gd.igc;
    gd.jstart = gd.jgc;
    gd.kstart = gd.kgc;

    gd.iend = gd.imax + gd.igc;
    gd.jend = gd.jmax + gd.jgc;
    gd.kend = gd.kmax + gd.kgc;

    check_ghost_cells();

    // allocate all arrays
    gd.x    .resize(gd.imax+2*gd.igc);
    gd.xh   .resize(gd.imax+2*gd.igc);
    gd.y    .resize(gd.jmax+2*gd.jgc);
    gd.yh   .resize(gd.jmax+2*gd.jgc);
    gd.z    .resize(gd.kmax+2*gd.kgc);
    gd.zh   .resize(gd.kmax+2*gd.kgc);
    gd.dz   .resize(gd.kmax+2*gd.kgc);
    gd.dzh  .resize(gd.kmax+2*gd.kgc);
    gd.dzi  .resize(gd.kmax+2*gd.kgc);
    gd.dzhi .resize(gd.kmax+2*gd.kgc);
    gd.dzi4 .resize(gd.kmax+2*gd.kgc);
    gd.dzhi4.resize(gd.kmax+2*gd.kgc);

    // initialize the communication functions
    init_mpi();

    // Initialize the transposes.
    transpose.init();
}

/**
 * This function initializes the fields containing the grid dimensions based
 * on the profiles in the input file.
 * @param inputin Pointer to the input class.
 */
template<typename TF>
void Grid<TF>::create(Data_block& profs)
{
    // Get the grid coordinates from the input.
    profs.get_vector(gd.z, "z", gd.kmax, 0, gd.kstart);

    if (gd.z[gd.kend-1] > gd.zsize)
    {
        master.print_error("Highest grid point is above prescribed zsize\n");
        throw 1;
    }

    // calculate the grid
    calculate();
}

/**
 * This function calculates the scalars and arrays that contain the information
 * on the grid spacing.
 */
template<typename TF>
void Grid<TF>::calculate()
{
    int i,j,k;

    // calculate the grid spacing
    gd.dx  = gd.xsize / gd.itot;
    gd.dy  = gd.ysize / gd.jtot;
    gd.dxi = 1./gd.dx;
    gd.dyi = 1./gd.dy;

    // calculate the offset per process to get the true x- and y-coordinate
    double xoff = master.mpicoordx * gd.xsize / master.npx;
    double yoff = master.mpicoordy * gd.ysize / master.npy;

    // calculate the x and y coordinates
    for (i=0; i<gd.icells; ++i)
    {
        gd.x [i] = 0.5*gd.dx + (i-gd.igc)*gd.dx + xoff;
        gd.xh[i] = (i-gd.igc)*gd.dx + xoff;
    }

    for (j=0; j<gd.jcells; ++j)
    {
        gd.y [j] = 0.5*gd.dy + (j-gd.jgc)*gd.dy + yoff;
        gd.yh[j] = (j-gd.jgc)*gd.dy + yoff;
    }

    // the calculation of ghost cells and flux levels has to go according to numerical scheme
    if (swspatialorder == "2")
    {
        gd.z[gd.kstart-1] = -gd.z[gd.kstart];
        gd.z[gd.kend]     = 2.*gd.zsize - gd.z[gd.kend-1];

        for (k=gd.kstart+1; k<gd.kend; ++k)
            gd.zh[k] = 0.5*(gd.z[k-1]+gd.z[k]);
        gd.zh[gd.kstart] = 0.;
        gd.zh[gd.kend]   = gd.zsize;

        // calculate the half levels according to the numerical scheme
        // compute the height of the grid cells
        for (k=1; k<gd.kcells; ++k)
        {
            gd.dzh [k] = gd.z[k] - gd.z[k-1];
            gd.dzhi[k] = 1./gd.dzh[k];
        }
        gd.dzh [gd.kstart-1] = gd.dzh [gd.kstart+1];
        gd.dzhi[gd.kstart-1] = gd.dzhi[gd.kstart+1];

        // compute the height of the grid cells
        for (k=1; k<gd.kcells-1; ++k)
        {
            gd.dz [k] = gd.zh[k+1] - gd.zh[k];
            gd.dzi[k] = 1./gd.dz[k];
        }
        gd.dz [gd.kstart-1] = gd.dz [gd.kstart];
        gd.dzi[gd.kstart-1] = gd.dzi[gd.kstart];
        gd.dz [gd.kend]     = gd.dz [gd.kend-1];
        gd.dzi[gd.kend]     = gd.dzi[gd.kend-1];

        // do not calculate 4th order gradients for 2nd order
    }

    if (swspatialorder == "4")
    {
        using namespace Finite_difference::O4;

        // calculate the height of the ghost cell
        gd.z[gd.kstart-1] = -2.*gd.z[gd.kstart] + (1./3.)*gd.z[gd.kstart+1];
        gd.z[gd.kstart-2] = -9.*gd.z[gd.kstart] +      2.*gd.z[gd.kstart+1];

        gd.z[gd.kend  ] = (8./3.)*gd.zsize - 2.*gd.z[gd.kend-1] + (1./3.)*gd.z[gd.kend-2];
        gd.z[gd.kend+1] =      8.*gd.zsize - 9.*gd.z[gd.kend-1] +      2.*gd.z[gd.kend-2];

        // Initialize the non-used values at a large value
        gd.z[gd.kstart-3] = Constants::dhuge;
        gd.z[gd.kend+2  ] = Constants::dhuge;

        gd.zh[gd.kstart] = 0.;
        for (k=gd.kstart+1; k<gd.kend; ++k)
            gd.zh[k] = ci0*gd.z[k-2] + ci1*gd.z[k-1] + ci2*gd.z[k] + ci3*gd.z[k+1];
        gd.zh[gd.kend] = gd.zsize;

        gd.zh[gd.kstart-1] = bi0*gd.z[gd.kstart-2] + bi1*gd.z[gd.kstart-1] + bi2*gd.z[gd.kstart] + bi3*gd.z[gd.kstart+1];
        gd.zh[gd.kend+1]   = ti0*gd.z[gd.kend-2  ] + ti1*gd.z[gd.kend-1  ] + ti2*gd.z[gd.kend  ] + ti3*gd.z[gd.kend+1  ];

        // calculate the half levels according to the numerical scheme
        // compute the height of the grid cells
        for (k=1; k<gd.kcells; ++k)
        {
            gd.dzh [k] = gd.z[k] - gd.z[k-1];
            gd.dzhi[k] = 1./gd.dzh[k];
        }
        gd.dzh [gd.kstart-3] = gd.dzh [gd.kstart+3];
        gd.dzhi[gd.kstart-3] = gd.dzhi[gd.kstart+3];

        // compute the height of the grid cells
        for (k=1; k<gd.kcells-1; ++k)
        {
            gd.dz [k] = gd.zh[k+1] - gd.zh[k];
            gd.dzi[k] = 1./gd.dz[k];
        }
        gd.dz [gd.kstart-3] = gd.dz [gd.kstart+2];
        gd.dzi[gd.kstart-3] = gd.dzi[gd.kstart+2];
        gd.dz [gd.kend+2] = gd.dz [gd.kend-3];
        gd.dzi[gd.kend+2] = gd.dzi[gd.kend-3];

        // calculate the fourth order gradients
        for (k=gd.kstart; k<gd.kend; ++k)
        {
            gd.dzi4 [k] = 1./(cg0*gd.zh[k-1] + cg1*gd.zh[k  ] + cg2*gd.zh[k+1] + cg3*gd.zh[k+2]);
            gd.dzhi4[k] = 1./(cg0*gd.z [k-2] + cg1*gd.z [k-1] + cg2*gd.z [k  ] + cg3*gd.z [k+1]);
        }
        gd.dzhi4[gd.kend] = 1./(cg0*gd.z[gd.kend-2] + cg1*gd.z[gd.kend-1] + cg2*gd.z[gd.kend] + cg3*gd.z[gd.kend+1]);

        // bc's
        gd.dzi4 [gd.kstart-1] = 1./(bg0*gd.zh[gd.kstart-1] + bg1*gd.zh[gd.kstart  ] + bg2*gd.zh[gd.kstart+1] + bg3*gd.zh[gd.kstart+2]);
        gd.dzhi4[gd.kstart-1] = 1./(bg0*gd.z [gd.kstart-2] + bg1*gd.z [gd.kstart-1] + bg2*gd.z [gd.kstart  ] + bg3*gd.z [gd.kstart+1]);

        gd.dzi4 [gd.kend  ] = 1./(tg0*gd.zh[gd.kend-2] + tg1*gd.zh[gd.kend-1] + tg2*gd.zh[gd.kend] + tg3*gd.zh[gd.kend+1]);
        gd.dzhi4[gd.kend+1] = 1./(tg0*gd.z [gd.kend-2] + tg1*gd.z [gd.kend-1] + tg2*gd.z [gd.kend] + tg3*gd.z [gd.kend+1]);

        // Define gradients at the boundary for the divgrad calculations.
        gd.dzhi4bot = 1./(bg0*gd.z[gd.kstart-1] + bg1*gd.z[gd.kstart] + bg2*gd.z[gd.kstart+1] + bg3*gd.z[gd.kstart+2]);
        gd.dzhi4top = 1./(tg0*gd.z[gd.kend-3  ] + tg1*gd.z[gd.kend-2] + tg2*gd.z[gd.kend-1  ] + tg3*gd.z[gd.kend    ]);

        // Initialize the unused values at a huge value to allow for easier error tracing.
        gd.dzi4[gd.kstart-2] = Constants::dhuge;
        gd.dzi4[gd.kstart-3] = Constants::dhuge;
        gd.dzi4[gd.kend+1  ] = Constants::dhuge;
        gd.dzi4[gd.kend+2  ] = Constants::dhuge;
    }
}

template<typename TF>
const Grid_data<TF>& Grid<TF>::get_grid_data()
{
    return gd;
}

template<typename TF>
void Grid<TF>::save()
{
    save_grid();
}

template<typename TF>
void Grid<TF>::load()
{
    load_grid();
}

/**
 * This function checks whether the number of ghost cells does not exceed the slice thickness.
 */
template<typename TF>
void Grid<TF>::check_ghost_cells()
{
    // Check whether the size per patch is larger than number of ghost cells for 3D runs.
    if (gd.imax < gd.igc)
    {
	    master.print_error("Patch size in x-dir (%d) is smaller than the number of ghost cells (%d).\n",(gd.iend-gd.istart), gd.igc);
	    master.print_error("Either increase itot or decrease npx.\n");
        throw 1;
    }

    // Check the jtot > 1 condition, to still allow for 2d runs.
    if (gd.jtot > 1 && gd.jmax < gd.jgc)
    {
	    master.print_error("Patch size in y-dir (%d) is smaller than the number of ghost cells (%d).\n",(gd.jend-gd.jstart), gd.jgc);
	    master.print_error("Either increase jtot or decrease npy.\n");
        throw 1;
    }
}

/**
 * This function increases the number of ghost cells in case necessary.
 * @param igc Ghost cells in the x-direction.
 * @param jgc Ghost cells in the y-direction.
 * @param kgc Ghost cells in the z-direction.
 */
template<typename TF>
void Grid<TF>::set_minimum_ghost_cells(const int igcin, const int jgcin, const int kgcin)
{
    gd.igc = std::max(gd.igc, igcin);
    gd.jgc = std::max(gd.jgc, jgcin);
    gd.kgc = std::max(gd.kgc, kgcin);

    // BvS: this doesn't work; imax is undefined if this routine is called from a class constructor
    // Removed it since this check is anyhow always performed from the init() of grid (after defining imax)
    //check_ghost_cells();
}

/**
 * This function does a second order horizontal interpolation in the x-direction
 * to the selected location on the grid.
 * @param out Pointer to the output field.
 * @param in Pointer to the input field.
 * @param locx Integer containing the location of the input field,
 * where a value of 1 refers to the flux level.
 */
template<typename TF>
void Grid<TF>::interpolate_2nd(TF* const restrict out, const TF* const restrict in, const int locin[3], const int locout[3])
{
    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const int iih = (locin[0]-locout[0])*ii;
    const int jjh = (locin[1]-locout[1])*jj;

    // interpolate the field
    // \TODO add the vertical component
    for (int k=0; k<gd.kcells; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                out[ijk] = 0.5*(0.5*in[ijk    ] + 0.5*in[ijk+iih    ])
                         + 0.5*(0.5*in[ijk+jjh] + 0.5*in[ijk+iih+jjh]);
            }
}

// /**
//  * This function does a fourth order horizontal interpolation in the x-direction
//  * to the selected location on the grid.
//  * @param out Pointer to the output field.
//  * @param in Pointer to the input field.
//  * @param locx Integer containing the location of the input field,
//  * where a value of 1 refers to the flux level.
//  */
// template<typename TF>
// void Grid<TF>::interpolate_4th(double* restrict out, double* restrict in, const int locin[3], const int locout[3])
// {
//     using namespace Finite_difference::O4;
//
//     // interpolation function, locx = 1 indicates that the reference is at the half level
//     const int ii = 1;
//     const int jj = icells;
//     const int kk = ijcells;
//
//     // a shift to the left gives minus 1 a shift to the right +1
//     const int iih1 = 1*(locin[0]-locout[0])*ii;
//     const int iih2 = 2*(locin[0]-locout[0])*ii;
//     const int jjh1 = 1*(locin[1]-locout[1])*jj;
//     const int jjh2 = 2*(locin[1]-locout[1])*jj;
//
//     // \TODO add the vertical component
//     for (int k=0; k<kcells; ++k)
//         for (int j=jstart; j<jend; ++j)
// #pragma ivdep
//             for (int i=istart; i<iend; ++i)
//             {
//                 const int ijk = i + j*jj + k*kk;
//                 out[ijk] = ci0*(ci0*in[ijk-iih1-jjh1] + ci1*in[ijk-jjh1] + ci2*in[ijk+iih1-jjh1] + ci3*in[ijk+iih2-jjh1])
//                          + ci1*(ci0*in[ijk-iih1     ] + ci1*in[ijk     ] + ci2*in[ijk+iih1     ] + ci3*in[ijk+iih2     ])
//                          + ci2*(ci0*in[ijk-iih1+jjh1] + ci1*in[ijk+jjh1] + ci2*in[ijk+iih1+jjh1] + ci3*in[ijk+iih2+jjh1])
//                          + ci3*(ci0*in[ijk-iih1+jjh2] + ci1*in[ijk+jjh2] + ci2*in[ijk+iih1+jjh2] + ci3*in[ijk+iih2+jjh2]);
//             }
// }
//
// template<typename TF>
// void Grid<TF>::calc_mean(double* restrict prof, const double* restrict data, const int krange)
// {
//     const int jj = icells;
//     const int kk = ijcells;
//
//     for (int k=0; k<krange; ++k)
//     {
//         prof[k] = 0.;
//         for (int j=jstart; j<jend; ++j)
// #pragma ivdep
//             for (int i=istart; i<iend; ++i)
//             {
//                 const int ijk  = i + j*jj + k*kk;
//                 prof[k] += data[ijk];
//             }
//     }
//
//     master.sum(prof, krange);
//
//     const double n = itot*jtot;
//
//     for (int k=0; k<krange; ++k)
//         prof[k] /= n;
// }

template class Grid<double>;
template class Grid<float>;
