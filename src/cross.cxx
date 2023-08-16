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
#include <cstdlib>
#include <cmath>
#include <algorithm>    // std::count
#include <stdexcept>

#include "master.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "cross.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "timeloop.h"

namespace
{
    template<typename TF>
    void calc_lngrad_4th(
            const TF* const restrict a, TF* const restrict lngrad,
            TF dxi, TF dyi, const TF* const restrict dzi4,
            int icells, int ijcells,
            int istart, int iend,
            int jstart, int jend,
            int kstart, int kend)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int jj3 = 3*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        // calculate the log of the gradient
        // bottom
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + kstart*kk1;
                lngrad[ijk] = std::log( Constants::dtiny +
                    std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-ii3] + ci1<TF>*a[ijk-ii2] + ci2<TF>*a[ijk-ii1] + ci3<TF>*a[ijk    ])
                              + cg1<TF>*(ci0<TF>*a[ijk-ii2] + ci1<TF>*a[ijk-ii1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+ii1])
                              + cg2<TF>*(ci0<TF>*a[ijk-ii1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+ii1] + ci3<TF>*a[ijk+ii2])
                              + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+ii1] + ci2<TF>*a[ijk+ii2] + ci3<TF>*a[ijk+ii3]) ) * dxi, 2.)

                  + std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-jj3] + ci1<TF>*a[ijk-jj2] + ci2<TF>*a[ijk-jj1] + ci3<TF>*a[ijk    ])
                              + cg1<TF>*(ci0<TF>*a[ijk-jj2] + ci1<TF>*a[ijk-jj1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+jj1])
                              + cg2<TF>*(ci0<TF>*a[ijk-jj1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+jj1] + ci3<TF>*a[ijk+jj2])
                              + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+jj1] + ci2<TF>*a[ijk+jj2] + ci3<TF>*a[ijk+jj3]) ) * dyi, 2.)

                  + std::pow( ( cg0<TF>*(bi0<TF>*a[ijk-kk2] + bi1<TF>*a[ijk-kk1] + bi2<TF>*a[ijk    ] + bi3<TF>*a[ijk+kk1])
                              + cg1<TF>*(ci0<TF>*a[ijk-kk2] + ci1<TF>*a[ijk-kk1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+kk1])
                              + cg2<TF>*(ci0<TF>*a[ijk-kk1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+kk1] + ci3<TF>*a[ijk+kk2])
                              + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+kk1] + ci2<TF>*a[ijk+kk2] + ci3<TF>*a[ijk+kk3]) ) * dzi4[kstart], 2.) );
            }

        // interior
        for (int k=kstart+1; k<kend-1; k++)
            for (int j=jstart; j<jend; j++)
    #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    lngrad[ijk] = std::log( Constants::dtiny +
                        std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-ii3] + ci1<TF>*a[ijk-ii2] + ci2<TF>*a[ijk-ii1] + ci3<TF>*a[ijk    ])
                                  + cg1<TF>*(ci0<TF>*a[ijk-ii2] + ci1<TF>*a[ijk-ii1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+ii1])
                                  + cg2<TF>*(ci0<TF>*a[ijk-ii1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+ii1] + ci3<TF>*a[ijk+ii2])
                                  + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+ii1] + ci2<TF>*a[ijk+ii2] + ci3<TF>*a[ijk+ii3]) ) * dxi, 2.)

                      + std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-jj3] + ci1<TF>*a[ijk-jj2] + ci2<TF>*a[ijk-jj1] + ci3<TF>*a[ijk    ])
                                  + cg1<TF>*(ci0<TF>*a[ijk-jj2] + ci1<TF>*a[ijk-jj1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+jj1])
                                  + cg2<TF>*(ci0<TF>*a[ijk-jj1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+jj1] + ci3<TF>*a[ijk+jj2])
                                  + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+jj1] + ci2<TF>*a[ijk+jj2] + ci3<TF>*a[ijk+jj3]) ) * dyi, 2.)

                      + std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-kk3] + ci1<TF>*a[ijk-kk2] + ci2<TF>*a[ijk-kk1] + ci3<TF>*a[ijk    ])
                                  + cg1<TF>*(ci0<TF>*a[ijk-kk2] + ci1<TF>*a[ijk-kk1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+kk1])
                                  + cg2<TF>*(ci0<TF>*a[ijk-kk1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+kk1] + ci3<TF>*a[ijk+kk2])
                                  + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+kk1] + ci2<TF>*a[ijk+kk2] + ci3<TF>*a[ijk+kk3]) ) * dzi4[k], 2.) );
                }

        // top
        for (int j=jstart; j<jend; j++)
        #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                lngrad[ijk] = std::log(Constants::dtiny +
                    std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-ii3] + ci1<TF>*a[ijk-ii2] + ci2<TF>*a[ijk-ii1] + ci3<TF>*a[ijk    ])
                              + cg1<TF>*(ci0<TF>*a[ijk-ii2] + ci1<TF>*a[ijk-ii1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+ii1])
                              + cg2<TF>*(ci0<TF>*a[ijk-ii1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+ii1] + ci3<TF>*a[ijk+ii2])
                              + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+ii1] + ci2<TF>*a[ijk+ii2] + ci3<TF>*a[ijk+ii3]) ) * dxi, 2.)

                  + std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-jj3] + ci1<TF>*a[ijk-jj2] + ci2<TF>*a[ijk-jj1] + ci3<TF>*a[ijk    ])
                              + cg1<TF>*(ci0<TF>*a[ijk-jj2] + ci1<TF>*a[ijk-jj1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+jj1])
                              + cg2<TF>*(ci0<TF>*a[ijk-jj1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+jj1] + ci3<TF>*a[ijk+jj2])
                              + cg3<TF>*(ci0<TF>*a[ijk    ] + ci1<TF>*a[ijk+jj1] + ci2<TF>*a[ijk+jj2] + ci3<TF>*a[ijk+jj3]) ) * dyi, 2.)

                  + std::pow( ( cg0<TF>*(ci0<TF>*a[ijk-kk3] + ci1<TF>*a[ijk-kk2] + ci2<TF>*a[ijk-kk1] + ci3<TF>*a[ijk    ])
                              + cg1<TF>*(ci0<TF>*a[ijk-kk2] + ci1<TF>*a[ijk-kk1] + ci2<TF>*a[ijk    ] + ci3<TF>*a[ijk+kk1])
                              + cg2<TF>*(ci0<TF>*a[ijk-kk1] + ci1<TF>*a[ijk    ] + ci2<TF>*a[ijk+kk1] + ci3<TF>*a[ijk+kk2])
                              + cg3<TF>*(ti0<TF>*a[ijk-kk1] + ti1<TF>*a[ijk    ] + ti2<TF>*a[ijk+kk1] + ti3<TF>*a[ijk+kk2]) ) * dzi4[kend-1], 2.) );
            }

    }

    template<typename TF>
    void calc_lngrad_2nd(
            const TF* const restrict a, TF* const restrict lngrad,
            TF dxi, TF dyi, const TF* const restrict dzi,
            int icells, int ijcells,
            int istart, int iend,
            int jstart, int jend,
            int kstart, int kend)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        // calculate the log of the gradient
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    lngrad[ijk] = std::log( Constants::dtiny
                            + std::pow( ( interp2(a[ijk   ], a[ijk+ii])
                                        - interp2(a[ijk-ii], a[ijk   ]) ) * dxi, TF(2))

                            + std::pow( ( interp2(a[ijk   ], a[ijk+jj])
                                        - interp2(a[ijk-jj], a[ijk   ]) ) * dyi, TF(2))

                            + std::pow( ( interp2(a[ijk   ], a[ijk+kk])
                                        - interp2(a[ijk-kk], a[ijk   ]) ) * dzi[k], TF(2)) );
                }
    }

    template<typename TF>
    void calc_cross_path(
            const TF* const restrict data, TF* const restrict tmp,
            const TF* const restrict rhoref, const TF* const restrict dz,
            int jj, int kk,
            int istart, int iend,
            int jstart, int jend,
            int kstart, int kend)
    {
    // Path is integrated in first full level, set to zero first
    for (int j=jstart; j<jend; j++)
        #pragma ivdep
        for (int i=istart; i<iend; i++)
        {
            const int ijk = i + j*jj + kstart*kk;
            tmp[ijk] = 0.;
        }

    // Integrate with height
    for (int k=kstart; k<kend; k++)
        for (int j=jstart; j<jend; j++)
        #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk1 = i + j*jj + kstart*kk;
                const int ijk  = i + j*jj + k*kk;
                tmp[ijk1] += rhoref[k] * data[ijk] * dz[k];
            }
    }

    template<typename TF>
    void calc_cross_height_threshold(
            const TF* const restrict data, TF* const restrict height,
            const TF* const restrict z, TF threshold, bool upward, TF fillvalue,
            int jj, int kk,
            int istart, int iend,
            int jstart, int jend,
            int kstart, int kend)
    {
        // Set height to NetCDF fill value
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jj;
                height[ij] = fillvalue;
            }

        if(upward) // Find lowest grid level where data > threshold
        {

            for (int j=jstart; j<jend; j++)
                for (int i=istart; i<iend; i++)
                    for (int k=kstart; k<kend; k++)
                    {
                        const int ij   = i + j*jj;
                        const int ijk  = i + j*jj + k*kk;

                        if(data[ijk] > threshold)
                        {
                            height[ij] = z[k];
                            break;
                        }
                    }
        }
        else // Find highest grid level where data > threshold
        {
            for (int j=jstart; j<jend; j++)
                for (int i=istart; i<iend; i++)
                    for (int k=kend-1; k>kstart-1; k--)
                    {
                        const int ij   = i + j*jj;
                        const int ijk  = i + j*jj + k*kk;

                        if(data[ijk] > threshold)
                        {
                            height[ij] = z[k];
                            break;
                        }
                    }
        }
    }
}

template<typename TF>
Cross<TF>::Cross(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), soil_grid(soilgridin),
    fields(fieldsin), field3d_io(master, grid)
{
    swcross = inputin.get_item<bool>("cross", "swcross", "", false);

    if (swcross)
    {
       // Get the time at which the cross sections are triggered.
        sampletime = inputin.get_item<double>("cross", "sampletime", "");

        // Get list of cross variables.
        crosslist = inputin.get_list<std::string>("cross", "crosslist", "", std::vector<std::string>());

        // Crash on empty list.
        if (crosslist.empty())
        {
            std::string msg = "Empty cross-section list";
            throw std::runtime_error(msg);
        }

        // Get the list of locations at which to take cross sections
        xy = inputin.get_list<TF>("cross", "xy", "", std::vector<TF>());
        xz = inputin.get_list<TF>("cross", "xz", "", std::vector<TF>());
        yz = inputin.get_list<TF>("cross", "yz", "", std::vector<TF>());

        // Get the list of vertical soil locations
        xy_soil = inputin.get_list<TF>("cross", "xy_soil", "", std::vector<TF>());
    }
}

template<typename TF>
Cross<TF>::~Cross()
{
}

// check whether saving the slice was successful and print appropriate message
template<typename TF>
int Cross<TF>::check_save(int error, char* filename)
{
    if (error == 0)
        return 0;
    else
    {
        master.print_message("Saving \"%s\" ... FAILED\n", filename);
        return 1;
    }
}

template<typename TF>
void Cross<TF>::init(double ifactor)
{
    if (!swcross)
        return;

    isampletime = static_cast<unsigned long>(ifactor * sampletime);
}

template<typename TF>
void Cross<TF>::create()
{
    int temploc, temploch, hoffset;

    auto& gd  = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();
    auto& md  = master.get_MPI_data();

    // Find nearest full and half grid locations of xz cross-sections.
    for (auto& it: xz)
    {
        // Find the index of the slice.
        temploc  = (int) floor(it/(gd.dy));
        temploch = (int) floor((it+(gd.dy/2.))/(gd.dy));

        if (it < 0 || it > gd.ysize) // Check if cross location is inside domain
        {
            std::string msg = std::to_string(it) + " in [cross][xz] is outside domain";
            throw std::runtime_error(msg);
        }
        else
        {
            if (it == gd.ysize) // Exception for full level when requesting domain size
                --temploc;

            // Find the corresponding index, make sure to handle MPI properly.
            TF ycross = -1.;
            if (temploc / gd.jmax == md.mpicoordy)
                ycross = gd.y[temploc % gd.jmax + gd.jgc];
            master.max(&ycross, 1);

            TF ycrossh = -1.;
            if (temploch / gd.jmax == md.mpicoordy)
                ycrossh = gd.y[temploch % gd.jmax + gd.jgc];
            master.max(&ycrossh, 1);

            if (std::find(jxz.begin(), jxz.end(), temploc) != jxz.end()) // Check for duplicate entries
                master.print_warning("Removed duplicate entry y=%f for [cross][xz]=%f\n", ycross, it);
            else // Add to cross-list
            {
                jxz.push_back(temploc);
                master.print_message("Added XZ cross at y=%f (j=%i) for [cross][xz]=%f\n", ycross, temploc, it);
            }

            if (std::find(jxzh.begin(), jxzh.end(), temploch) != jxzh.end()) // Check for duplicate entries
                master.print_warning("Removed duplicate entry yh=%f for [cross][xz]=%f\n", ycrossh, it);
            else // Add to cross-list
            {
                jxzh.push_back(temploch);
                master.print_message("Added XZ cross at yh=%f (j=%i) for [cross][xz]=%f\n", ycrossh, temploch, it);
            }
        }
    }

    // Find nearest full and half grid locations of yz cross-sections.
    for (auto& it: yz)
    {
        temploc  = (int) floor(it/(gd.dx));
        temploch = (int) floor((it+(gd.dx/2.))/(gd.dx));

        if (it < 0 || it > gd.xsize) // Check if cross location is inside domain
        {
            std::string msg = std::to_string(it) + " in [cross][yz] is outside domain";
            throw std::runtime_error(msg);
        }
        else
        {
            if (it == gd.xsize) // Exception for full level when requesting domain size
                --temploc;

            // Find the corresponding index, make sure to handle MPI properly.
            TF xcross = -1.;
            if (temploc / gd.imax == md.mpicoordx)
                xcross = gd.x[temploc % gd.imax + gd.igc];
            master.max(&xcross, 1);

            TF xcrossh = -1.;
            if (temploch / gd.imax == md.mpicoordx)
                xcrossh = gd.x[temploch % gd.imax + gd.igc];
            master.max(&xcrossh, 1);

            if (std::find(ixz.begin(), ixz.end(), temploc) != ixz.end()) // Check for duplicate entries
                master.print_warning("Removed duplicate entry x=%f for [cross][yz]=%f\n", xcross, it);
            else // Add to cross-list
            {
                ixz.push_back(temploc);
                master.print_message("Added YZ cross at x=%f (i=%i) for [cross][yz]=%f\n", xcross, temploc, it);
            }

            if (std::find(ixzh.begin(), ixzh.end(), temploch) != ixzh.end()) // Check for duplicate entries
                master.print_warning("Removed duplicate entry xh=%f for [cross][yz]=%f\n", xcrossh, it);
            else // Add to cross-list
            {
                ixzh.push_back(temploch);
                master.print_message("Added YZ cross at xh=%f (i=%i) for [cross][yz]=%f\n", xcrossh, temploch, it);
            }
        }
    }

    // Find nearest full and half grid locations of xy cross-sections.
    for (auto& it: xy)
    {
        hoffset = 0;
        if (it < 0 || it > gd.zsize) // Check if cross location is inside domain
        {
            std::string msg = std::to_string(it) + " in [cross][xy] is outside domain";
            throw std::runtime_error(msg);
        }
        else
        {
            if (it == gd.zsize) // Exception for domain top: use half level at domain top, full level below
            {
                temploc = gd.kmax-1;
                hoffset = 1;
            }
            else
            {
                for (int k=gd.kstart; k<gd.kend; k++) // Loop over height to find the nearest full level
                {
                    if ((it >= gd.zh[k]) && (it < gd.zh[k+1]))
                    {
                        temploc = k - gd.kgc;
                        if (it >= gd.z[k]) // Add offset for half level
                            hoffset = 1;
                        break;
                    }
                }
            }

            if (std::find(kxy.begin(), kxy.end(), temploc) != kxy.end()) // Check for duplicate entries
                master.print_warning("Removed duplicate entry z=%f for [cross][xy]=%f\n", gd.z[temploc+gd.kgc],it);
            else // Add to cross-list
            {
                kxy.push_back(temploc);
                master.print_message("Added XY cross at z=%f (k=%i) for [cross][xy]=%f\n", gd.z[temploc+gd.kgc],temploc,it);
            }

            if (std::find(kxyh.begin(), kxyh.end(), temploc+hoffset) != kxyh.end()) // Check for duplicate entries
                master.print_warning("Removed duplicate entry zh=%f for [cross][xy]=%f\n", gd.zh[temploc+hoffset+gd.kgc],it);
            else // Add to cross-list
            {
                kxyh.push_back(temploc+hoffset);
                master.print_message("Added XY cross at zh=%f (k=%i) for [cross][xy]=%f\n", gd.zh[temploc+hoffset+gd.kgc],temploc+hoffset,it);
            }
        }
    }

    // Find nearest index of XY cross-sections in soil
    for (auto& z : xy_soil)
    {
        if (z > 0 || z < sgd.zsize)
        {
            std::string msg = std::to_string(z) + " in [cross][xy_soil] is outside soil";
            throw std::runtime_error(msg);
        }
        else
        {
            for (int k=sgd.kstart; k<sgd.kend; ++k)
            {
                if ((z >= sgd.zh[k]) && (z < sgd.zh[k+1]))
                {
                    temploc = k - sgd.kgc;
                    break;
                }
            }

            if (std::find(kxy_soil.begin(), kxy_soil.end(), temploc) != kxy_soil.end()) // Check for duplicate entries
                master.print_warning("Removed duplicate entry z=%f for [cross][xy_soil]=%f\n", sgd.z[temploc+sgd.kgc], z);
            else // Add to cross-list
            {
                kxy_soil.push_back(temploc);
                master.print_message("Added XY cross at z=%f (k=%i) for [cross][xy_soil]=%f\n", sgd.z[temploc+sgd.kgc],temploc, z);
            }
        }
    }

    /* All classes (fields, thermo, boundary) have removed their cross-variables from
       crosslist by now. If it isnt empty, print warnings for invalid variables */
    if (crosslist.size() > 0)
    {
    for (auto& it: crosslist)
            master.print_warning("field %s in [cross][crosslist] is illegal\n", it.c_str());
    }

}

template<typename TF>
unsigned long Cross<TF>::get_time_limit(unsigned long itime)
{
    if (!swcross)
        return Constants::ulhuge;

    unsigned long idtlim = isampletime - itime % isampletime;

    return idtlim;
}

template<typename TF>
bool Cross<TF>::do_cross(unsigned long itime)
{
    // check if cross are enabled
    if (!swcross)
        return false;


    // check if time for execution
    if (itime % isampletime != 0)
        return false;

    // return true such that cross are computed
    return true;
}


template<typename TF>
std::vector<std::string>& Cross<TF>::get_crosslist()
{
    return crosslist;
}

/*
This function is called from various classes, "allowed_variables" is a vector
with the variables that each class can potentially provide, get_enabled_variables()
cross-references this with the cross-sections requested in the .ini file (crosslist),
removes the cross-sections that the class can and should provide from crosslist,
and returns a vector containing those variables (difficult explanation...)
*/
template<typename TF>
std::vector<std::string> Cross<TF>::get_enabled_variables(const std::vector<std::string>& allowed_variables)
{
    std::vector<std::string> variables_to_return;

    std::vector<std::string>::iterator it=crosslist.begin();
    while (it != crosslist.end())
    {
        if (std::count(allowed_variables.begin(), allowed_variables.end(), *it))
        {
            // Remove variable from global list, put in list which is returned
            variables_to_return.push_back(*it);
            crosslist.erase(it); // erase() returns iterator of next element..
        }
        else
            ++it;
    }

    return variables_to_return;
}


template<typename TF>
int Cross<TF>::cross_simple(
        TF* restrict data, TF restrict offset, const std::string& name, const int iotime, const std::array<int,3>& loc)
{
    auto& gd = grid.get_grid_data();

    int nerror = 0;
    char filename[256];

    auto tmpfld = fields.get_tmp();
    auto tmp = tmpfld->fld.data();

    // Loop over the index arrays to save all xz cross sections.
    if (loc == gd.vloc)
    {
        for (auto& it: jxzh)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", it, iotime);
            nerror += check_save(
                    field3d_io.save_xz_slice(data, offset, tmp, filename, it, gd.kstart, gd.kend), filename);
        }
    }
    else
    {
        for (auto& it: jxz)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", it, iotime);
            nerror += check_save(
                    field3d_io.save_xz_slice(data, offset, tmp, filename, it, gd.kstart, gd.kend), filename);
        }
    }

    // Loop over the index arrays to save all yz cross sections.
    if (loc == gd.uloc)
    {
        for (auto& it: ixzh)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "yz", it, iotime);
            nerror += check_save(
                    field3d_io.save_yz_slice(data, offset, tmp, filename, it, gd.kstart, gd.kend), filename);
        }
    }
    else
    {
        for (auto& it: ixz)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "yz", it, iotime);
            nerror += check_save(
                    field3d_io.save_yz_slice(data, offset, tmp, filename, it, gd.kstart, gd.kend), filename);
        }
    }

    if (loc == gd.wloc)
    {
        // loop over the index arrays to save all xy cross sections
        for (auto& it: kxyh)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", it, iotime);
            nerror += check_save(field3d_io.save_xy_slice(data, offset, tmp, filename, it+gd.kgc), filename);
        }
    }
    else
    {
        for (auto& it: kxy)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", it, iotime);
            nerror += check_save(field3d_io.save_xy_slice(data, offset, tmp, filename, it+gd.kgc), filename);
        }
    }
    fields.release_tmp(tmpfld);

    return nerror;
}

template<typename TF>
int Cross<TF>::cross_plane(TF* restrict data, TF restrict offset, std::string name, int iotime)
{
    int nerror = 0;
    char filename[256];

    auto tmpfld = fields.get_tmp();
    auto tmp = tmpfld->fld.data();

    std::sprintf(filename, "%s.%s.%07d", name.c_str(), "xy", iotime);
    nerror += check_save(field3d_io.save_xy_slice(data, offset, tmp, filename), filename);
    fields.release_tmp(tmpfld);
    return nerror;
}

template<typename TF>
int Cross<TF>::cross_lngrad(TF* restrict a, std::string name, int iotime)
{
    auto& gd = grid.get_grid_data();

    int nerror = 0;
    char filename[256];

    auto lngradfld = fields.get_tmp();
    auto lngrad = lngradfld->fld.data();
    auto tmpfld = fields.get_tmp();
    auto tmp = tmpfld->fld.data();

    if (grid.get_spatial_order() == Grid_order::Second)
        calc_lngrad_2nd<TF>(
                a, lngrad, gd.dxi, gd.dyi, gd.dzi.data(),
                gd.icells, gd.ijcells, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend);
    else if (grid.get_spatial_order() == Grid_order::Fourth)
        calc_lngrad_4th<TF>(
                a, lngrad, gd.dxi, gd.dyi, gd.dzi4.data(),
                gd.icells, gd.ijcells, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend);

    // loop over the index arrays to save all xz cross sections
    TF no_offset = 0;
    for (auto& it: jxz)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", it, iotime);
        nerror += check_save(
                field3d_io.save_xz_slice(lngrad, no_offset, tmp, filename, it, gd.kstart, gd.kend),filename);
    }

    // loop over the index arrays to save all yz cross sections
    for (auto& it: ixz)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "yz", it, iotime);
        nerror += check_save(
                field3d_io.save_yz_slice(lngrad, no_offset, tmp, filename, it, gd.kstart, gd.kend),filename);
    }

    // loop over the index arrays to save all xy cross sections
    for (auto& it: kxy)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", it, iotime);
        nerror += check_save(field3d_io.save_xy_slice(lngrad, no_offset, tmp, filename, it+gd.kgc),filename);
    }

    fields.release_tmp(tmpfld);
    fields.release_tmp(lngradfld);

    return nerror;
}

template<typename TF>
int Cross<TF>::cross_path(TF* restrict data, std::string name, int iotime)
{

    int nerror = 0;
    TF no_offset = 0.;
    auto tmpfld = fields.get_tmp();
    auto tmp = tmpfld->fld.data();
    auto& gd = grid.get_grid_data();

    calc_cross_path<TF>(
            data, tmp, fields.rhoref.data(), gd.dz.data(),
            gd.icells, gd.ijcells,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend);

    nerror += cross_plane(&tmp[gd.kstart*gd.ijcells], no_offset, name, iotime);
    fields.release_tmp(tmpfld);
    return nerror;
}

/**
 * This routine calculates the lowest or highest height where data > threshold,
 * and writes a cross-section of the resulting height field
 * @param data Pointer to input data
 * @param height Pointer to 2D temporary field to store the height
 * @param tmp1 Pointer to temporary field for writing the cross-section
 * @param z Pointer to 1D field containing the levels of data
 * @param threshold Threshold value
 * @param Cross_direction Switch for bottom-up (Bottom_to_top) or top-down (Top_to_bottom)
 * @param name String containing the output name of the cross-section
 */
template<typename TF>
int Cross<TF>::cross_height_threshold(TF* restrict data, TF threshold, Cross_direction direction, std::string name, int iotime)
{

    auto& gd = grid.get_grid_data();
    int nerror = 0;
    TF no_offset = 0.;
    auto tmpfld = fields.get_tmp();
    auto height = tmpfld->fld.data();

    TF fillvalue = -1e-9; //TODO: SET FILL VALUE
    bool isupward = (direction == Cross_direction::Bottom_to_top);
    calc_cross_height_threshold<TF>(
            data, height, gd.z.data(), threshold, isupward, fillvalue,
            gd.icells, gd.ijcells,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend);

    nerror += cross_plane(height, no_offset, name, iotime);
    fields.release_tmp(tmpfld);
    return nerror;
}

template<typename TF>
int Cross<TF>::cross_soil(
        TF* restrict data, const std::string& name, const int iotime)
{
    auto& sgd = soil_grid.get_grid_data();

    int nerror = 0;
    char filename[256];
    TF no_offset = 0.;

    auto tmpfld = fields.get_tmp();
    auto tmp = tmpfld->fld.data();

    for (auto& it: jxz)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", it, iotime);
        nerror += check_save(
                field3d_io.save_xz_slice(data, no_offset, tmp, filename, it, sgd.kstart, sgd.kend), filename);
    }

    for (auto& it: ixz)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "yz", it, iotime);
        nerror += check_save(
                field3d_io.save_yz_slice(data, no_offset, tmp, filename, it, sgd.kstart, sgd.kend), filename);
    }

    for (auto& it: kxy_soil)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", it, iotime);
        nerror += check_save(field3d_io.save_xy_slice(data, no_offset, tmp, filename, it+sgd.kgc), filename);
    }

    fields.release_tmp(tmpfld);

    return nerror;
}


template class Cross<double>;
template class Cross<float>;
