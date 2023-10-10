/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
#include <iostream>
#include <algorithm>

#include "master.h"
#include "grid.h"
#include "input.h"
#include "netcdf_interface.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "timedep.h"
#include "stats.h"



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

    gd.utrans = input.get_item<TF>("grid", "utrans", "", 0.);
    gd.vtrans = input.get_item<TF>("grid", "vtrans", "", 0.);

    gd.lat = input.get_item<TF>("grid", "lat", "",  0.); 
    gd.lon = input.get_item<TF>("grid", "lon", "",  0.); 

    std::string swspatialorder = input.get_item<std::string>("grid", "swspatialorder", "");

    if (swspatialorder == "2")
        spatial_order = Grid_order::Second;
    else if (swspatialorder == "4")
        spatial_order = Grid_order::Fourth;
    else
    {
            std::string msg = swspatialorder + " is an illegal value for swspatialorder";
            throw std::runtime_error(msg);
    }

    // 2nd order scheme requires only 1 ghost cell
    if (spatial_order == Grid_order::Second)
    {
        gd.igc = 1;
        gd.jgc = 1;
        gd.kgc = 1;
    }
    // 4th order scheme requires 3 ghost cells
    else if (spatial_order == Grid_order::Fourth)
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
    auto& md = master.get_MPI_data();

    // Check whether the grid fits the processor configuration.
    if (gd.itot % md.npx != 0)
    {
        std::string msg = "itot = " + std::to_string(gd.itot) +  " is not a multiple of npx = " + std::to_string(md.npx);
        throw std::runtime_error(msg);
    }
    if (gd.itot % md.npy != 0)
    {
        std::string msg = "itot = " + std::to_string(gd.itot) +  " is not a multiple of npy = " + std::to_string(md.npy);
        throw std::runtime_error(msg);
    }
    // Check this one only when npy > 1, since the transpose in that direction only happens then.
    if (gd.jtot % md.npx != 0 && md.npy > 1)
    {
        std::string msg = "jtot = " + std::to_string(gd.jtot) +  " is not a multiple of npx = " + std::to_string(md.npx);
        throw std::runtime_error(msg);
    }
    if (gd.jtot % md.npy != 0)
    {
        std::string msg = "jtot = " + std::to_string(gd.jtot) +  " is not a multiple of npy = " + std::to_string(md.npy);
        throw std::runtime_error(msg);
    }
    if (gd.ktot % md.npx != 0)
    {
        std::string msg = "ktot = " + std::to_string(gd.ktot) +  " is not a multiple of npx = " + std::to_string(md.npx);
        throw std::runtime_error(msg);
    }

    // Calculate the total number of grid cells.
    gd.ntot = gd.itot*gd.jtot*gd.ktot;

    // Calculate the grid dimensions per process.
    gd.imax = gd.itot / md.npx;
    gd.jmax = gd.jtot / md.npy;
    gd.kmax = gd.ktot;
    gd.nmax = gd.imax*gd.jmax*gd.kmax;

    // Calculate the block sizes for the transposes.
    gd.iblock = gd.itot / md.npy;
    gd.jblock = gd.jtot / md.npx;
    gd.kblock = gd.ktot / md.npx;

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
void Grid<TF>::create(Input& inputin, Netcdf_handle& input_nc)
{
    // Get the grid coordinates from the input. This one is read from the global scope.
    input_nc.get_variable(gd.z, "z", {0}, {gd.ktot});
    std::rotate(gd.z.rbegin(), gd.z.rbegin() + gd.kstart, gd.z.rend());

    if (gd.z[gd.kend-1] > gd.zsize)
    {
        std::string msg = "Highest grid point is above prescribed zsize";
        throw std::runtime_error(msg);
    }
    // calculate the grid
    calculate();
}

template<typename TF>
void Grid<TF>::create_stats(Stats<TF>& stats)
{
    const std::string group_name = "default";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        stats.add_time_series("lat", "Latitude", "degrees", group_name);
        stats.add_time_series("lon", "Longitude", "degrees", group_name);
    }
}


template<typename TF>
void Grid<TF>::exec_stats(Stats<TF>& stats)
{
    stats.set_time_series("lat", gd.lat);
    stats.set_time_series("lon", gd.lon);
}
/**
 * This function calculates the scalars and arrays that contain the information
 * on the grid spacing.
 */
template<typename TF>
void Grid<TF>::calculate()
{
    auto& md = master.get_MPI_data();

    // calculate the grid spacing
    gd.dx  = gd.xsize / gd.itot;
    gd.dy  = gd.ysize / gd.jtot;
    gd.dxi = 1./gd.dx;
    gd.dyi = 1./gd.dy;

    // calculate the offset per process to get the true x- and y-coordinate
    double xoff = md.mpicoordx * gd.xsize / md.npx;
    double yoff = md.mpicoordy * gd.ysize / md.npy;

    // calculate the x and y coordinates
    for (int i=0; i<gd.icells; ++i)
    {
        gd.x [i] = 0.5*gd.dx + (i-gd.igc)*gd.dx + xoff;
        gd.xh[i] = (i-gd.igc)*gd.dx + xoff;
    }

    for (int j=0; j<gd.jcells; ++j)
    {
        gd.y [j] = 0.5*gd.dy + (j-gd.jgc)*gd.dy + yoff;
        gd.yh[j] = (j-gd.jgc)*gd.dy + yoff;
    }

    // the calculation of ghost cells and flux levels has to go according to numerical scheme
    if (spatial_order == Grid_order::Second)
    {
        gd.z[gd.kstart-1] = -gd.z[gd.kstart];
        gd.z[gd.kend]     = 2.*gd.zsize - gd.z[gd.kend-1];

        for (int k=gd.kstart+1; k<gd.kend; ++k)
            gd.zh[k] = 0.5*(gd.z[k-1]+gd.z[k]);
        gd.zh[gd.kstart] = 0.;
        gd.zh[gd.kend]   = gd.zsize;

        // calculate the half levels according to the numerical scheme
        // compute the height of the grid cells
        for (int k=1; k<gd.kcells; ++k)
        {
            gd.dzh [k] = gd.z[k] - gd.z[k-1];
            gd.dzhi[k] = 1./gd.dzh[k];
        }
        gd.dzh [gd.kstart-1] = gd.dzh [gd.kstart+1];
        gd.dzhi[gd.kstart-1] = gd.dzhi[gd.kstart+1];

        // compute the height of the grid cells
        for (int k=1; k<gd.kcells-1; ++k)
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

    if (spatial_order == Grid_order::Fourth)
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
        for (int k=gd.kstart+1; k<gd.kend; ++k)
            gd.zh[k] = ci0<TF>*gd.z[k-2] + ci1<TF>*gd.z[k-1] + ci2<TF>*gd.z[k] + ci3<TF>*gd.z[k+1];
        gd.zh[gd.kend] = gd.zsize;

        gd.zh[gd.kstart-1] = bi0<TF>*gd.z[gd.kstart-2] + bi1<TF>*gd.z[gd.kstart-1] + bi2<TF>*gd.z[gd.kstart] + bi3<TF>*gd.z[gd.kstart+1];
        gd.zh[gd.kend+1]   = ti0<TF>*gd.z[gd.kend-2  ] + ti1<TF>*gd.z[gd.kend-1  ] + ti2<TF>*gd.z[gd.kend  ] + ti3<TF>*gd.z[gd.kend+1  ];

        // calculate the half levels according to the numerical scheme
        // compute the height of the grid cells
        for (int k=1; k<gd.kcells; ++k)
        {
            gd.dzh [k] = gd.z[k] - gd.z[k-1];
            gd.dzhi[k] = 1./gd.dzh[k];
        }
        gd.dzh [gd.kstart-3] = gd.dzh [gd.kstart+3];
        gd.dzhi[gd.kstart-3] = gd.dzhi[gd.kstart+3];

        // compute the height of the grid cells
        for (int k=1; k<gd.kcells-1; ++k)
        {
            gd.dz [k] = gd.zh[k+1] - gd.zh[k];
            gd.dzi[k] = 1./gd.dz[k];
        }
        gd.dz [gd.kstart-3] = gd.dz [gd.kstart+2];
        gd.dzi[gd.kstart-3] = gd.dzi[gd.kstart+2];
        gd.dz [gd.kend+2] = gd.dz [gd.kend-3];
        gd.dzi[gd.kend+2] = gd.dzi[gd.kend-3];

        // calculate the fourth order gradients
        for (int k=gd.kstart; k<gd.kend; ++k)
        {
            gd.dzi4 [k] = 1./(cg0<TF>*gd.zh[k-1] + cg1<TF>*gd.zh[k  ] + cg2<TF>*gd.zh[k+1] + cg3<TF>*gd.zh[k+2]);
            gd.dzhi4[k] = 1./(cg0<TF>*gd.z [k-2] + cg1<TF>*gd.z [k-1] + cg2<TF>*gd.z [k  ] + cg3<TF>*gd.z [k+1]);
        }
        gd.dzhi4[gd.kend] = 1./(cg0<TF>*gd.z[gd.kend-2] + cg1<TF>*gd.z[gd.kend-1] + cg2<TF>*gd.z[gd.kend] + cg3<TF>*gd.z[gd.kend+1]);

        // bc's
        gd.dzi4 [gd.kstart-1] = 1./(bg0<TF>*gd.zh[gd.kstart-1] + bg1<TF>*gd.zh[gd.kstart  ] + bg2<TF>*gd.zh[gd.kstart+1] + bg3<TF>*gd.zh[gd.kstart+2]);
        gd.dzhi4[gd.kstart-1] = 1./(bg0<TF>*gd.z [gd.kstart-2] + bg1<TF>*gd.z [gd.kstart-1] + bg2<TF>*gd.z [gd.kstart  ] + bg3<TF>*gd.z [gd.kstart+1]);

        gd.dzi4 [gd.kend  ] = 1./(tg0<TF>*gd.zh[gd.kend-2] + tg1<TF>*gd.zh[gd.kend-1] + tg2<TF>*gd.zh[gd.kend] + tg3<TF>*gd.zh[gd.kend+1]);
        gd.dzhi4[gd.kend+1] = 1./(tg0<TF>*gd.z [gd.kend-2] + tg1<TF>*gd.z [gd.kend-1] + tg2<TF>*gd.z [gd.kend] + tg3<TF>*gd.z [gd.kend+1]);

        // Define gradients at the boundary for the divgrad calculations.
        gd.dzhi4bot = 1./(bg0<TF>*gd.z[gd.kstart-1] + bg1<TF>*gd.z[gd.kstart] + bg2<TF>*gd.z[gd.kstart+1] + bg3<TF>*gd.z[gd.kstart+2]);
        gd.dzhi4top = 1./(tg0<TF>*gd.z[gd.kend-3  ] + tg1<TF>*gd.z[gd.kend-2] + tg2<TF>*gd.z[gd.kend-1  ] + tg3<TF>*gd.z[gd.kend    ]);

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
void Grid<TF>::load(Input& inputin, Netcdf_handle& input_nc)
{
    load_grid();

    std::string timedep_dim = "time_latlon";
    swtimedep = inputin.get_item<bool>("grid", "swtimedep", "", false);
    tdep_latlon.emplace("lat", new Timedep<TF>(master, (*this), "lat", swtimedep));
    tdep_latlon.at("lat")->create_timedep(input_nc, timedep_dim);
    tdep_latlon.emplace("lon", new Timedep<TF>(master, (*this), "lon", swtimedep));
    tdep_latlon.at("lon")->create_timedep(input_nc, timedep_dim);


}

template <typename TF>
void Grid<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
        tdep_latlon.at("lat")->update_time_dependent(gd.lat, timeloop);
        tdep_latlon.at("lon")->update_time_dependent(gd.lon, timeloop);
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
        master.print_message("Patch size in x-dir (%d) is smaller than the number of ghost cells (%d).\n",(gd.iend-gd.istart), gd.igc);
        std::string msg = "Either increase itot or decrease npx";
        throw std::runtime_error(msg);
    }

    // Check the jtot > 1 condition, to still allow for 2d runs.
    if (gd.jtot > 1 && gd.jmax < gd.jgc)
    {
        master.print_message("Patch size in y-dir (%d) is smaller than the number of ghost cells (%d).\n",(gd.jend-gd.jstart), gd.jgc);
        std::string msg = "Either increase jtot or decrease npy";
        throw std::runtime_error(msg);
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
    const int kkh = (locin[2]-locout[2])*kk;

    // interpolate the field
    // \TODO add the vertical component
    int kstart = 0;
    int kend = gd.kcells;

    if (kkh == -kk)
        kstart = 1;
    if (kkh == kk)
        kend = gd.kcells - 1;

    #pragma omp parallel for
    for (int k=kstart; k<kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                out[ijk] = 0.25*(0.5*in[ijk            ] + 0.5*in[ijk+iih        ])
                         + 0.25*(0.5*in[ijk    +jjh    ] + 0.5*in[ijk+iih+jjh    ])
                         + 0.25*(0.5*in[ijk        +kkh] + 0.5*in[ijk+iih    +kkh])
                         + 0.25*(0.5*in[ijk    +jjh+kkh] + 0.5*in[ijk+iih+jjh+kkh]);
            }

    if (kkh == -kk) //ASSUMING ZERO AT BOTTOM BC
    {
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj;
                out[ijk] = 0.25*(0.5*in[ijk            ] + 0.5*in[ijk+iih        ])
                         + 0.25*(0.5*in[ijk    +jjh    ] + 0.5*in[ijk+iih+jjh    ]);
            }

    }
    else if (kkh == kk) //ASSUMING CONSTANT GRADIENT AT TOP BC
    {
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + kend*kk;
                out[ijk] = 1. *(0.5*in[ijk            ] + 0.5*in[ijk+iih        ])
                         + 1. *(0.5*in[ijk    +jjh    ] + 0.5*in[ijk+iih+jjh    ])
                         - 0.5*(0.5*in[ijk        -kkh] + 0.5*in[ijk+iih    -kkh])
                         - 0.5*(0.5*in[ijk    +jjh-kkh] + 0.5*in[ijk+iih+jjh-kkh]);
            }
    }
}

/**
 * This function does a fourth order horizontal interpolation in the x-direction
 * to the selected location on the grid.
 * @param out Pointer to the output field.
 * @param in Pointer to the input field.
 * @param locx Integer containing the location of the input field,
 * where a value of 1 refers to the flux level.
 */
template<typename TF>
void Grid<TF>::interpolate_4th(TF* restrict out, const TF* restrict in, const int locin[3], const int locout[3])
{
    using namespace Finite_difference::O4;

    // interpolation function, locx = 1 indicates that the reference is at the half level
    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    // a shift to the left gives minus 1 a shift to the right +1
    const int iih1 = 1*(locin[0]-locout[0])*ii;
    const int iih2 = 2*(locin[0]-locout[0])*ii;
    const int jjh1 = 1*(locin[1]-locout[1])*jj;
    const int jjh2 = 2*(locin[1]-locout[1])*jj;

    // \TODO add the vertical component
    #pragma omp parallel for
    for (int k=0; k<gd.kcells; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                out[ijk] = ci0<TF>*(ci0<TF>*in[ijk-iih1-jjh1] + ci1<TF>*in[ijk-jjh1] + ci2<TF>*in[ijk+iih1-jjh1] + ci3<TF>*in[ijk+iih2-jjh1])
                         + ci1<TF>*(ci0<TF>*in[ijk-iih1     ] + ci1<TF>*in[ijk     ] + ci2<TF>*in[ijk+iih1     ] + ci3<TF>*in[ijk+iih2     ])
                         + ci2<TF>*(ci0<TF>*in[ijk-iih1+jjh1] + ci1<TF>*in[ijk+jjh1] + ci2<TF>*in[ijk+iih1+jjh1] + ci3<TF>*in[ijk+iih2+jjh1])
                         + ci3<TF>*(ci0<TF>*in[ijk-iih1+jjh2] + ci1<TF>*in[ijk+jjh2] + ci2<TF>*in[ijk+iih1+jjh2] + ci3<TF>*in[ijk+iih2+jjh2]);
            }
}

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


#ifdef FLOAT_SINGLE
template class Grid<float>;
#else
template class Grid<double>;
#endif
