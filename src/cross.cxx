/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "cross.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"
#include "thermo.h"
#include "timeloop.h"
#include <netcdf>

Cross::Cross(Model* modelin, Input* inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    // Optional, by default switch cross off.
    int nerror = 0;
    nerror += inputin->get_item(&swcross, "cross", "swcross", "", "0");

    if (swcross == "1")
    {
        // Get the time at which the cross sections are triggered.
        nerror += inputin->get_item(&sampletime, "cross", "sampletime", "");

        // Get list of cross variables.
        nerror += inputin->get_list(&crosslist , "cross", "crosslist" , "");

        // Crash on empty list.
        if (crosslist.empty())
        {
            master->print_error("empty cross section list\n");
            throw 1;
        }

        // get the list of indices at which to take cross sections
        nerror += inputin->get_list(&xz, "cross", "xz", "");
        nerror += inputin->get_list(&yz, "cross", "yz", "");
        nerror += inputin->get_list(&xy, "cross", "xy", "");
    }

    if (nerror)
        throw 1;
}

Cross::~Cross()
{
}

// check whether saving the slice was successful and print appropriate message
int Cross::check_save(int error, char* filename)
{
    master->print_message("Saving \"%s\" ... ", filename);
    if (error == 0)
    {
        master->print_message("OK\n");
        return 0;
    }
    else
    {
        master->print_message("FAILED\n");
        return 1;
    }
}

void Cross::init(double ifactor)
{
    if (swcross == "0")
        return;

    isampletime = (unsigned long)(ifactor * sampletime);
}

void Cross::create()
{  
    int nerror = 0;
    int temploc, temploch, hoffset;

    // Find nearest full and half grid locations of xz cross-sections.
    for (std::vector<double>::iterator it=xz.begin(); it<xz.end(); ++it)
    {
        // Find the index of the slice.
        temploc  = (int) floor(*it/(grid->dy));
        temploch = (int) floor((*it+(grid->dy/2.))/(grid->dy));

        if (*it < 0 || *it > grid->ysize) // Check if cross location is inside domain
        {
            master->print_error("%f in [cross][xz] is outside domain\n", *it);
            ++nerror;
        }
        else
        {
            if (*it == grid->ysize) // Exception for full level when requesting domain size
                --temploc;

            // Find the corresponding index, make sure to handle MPI properly.
            double ycross = -1.;
            if (temploc / grid->jmax == master->mpicoordy)
                ycross = grid->y[temploc % grid->jmax + grid->jgc];
            master->max(&ycross, 1);

            double ycrossh = -1.;
            if (temploch / grid->jmax == master->mpicoordy)
                ycrossh = grid->y[temploch % grid->jmax + grid->jgc];
            master->max(&ycrossh, 1);

            if (std::find(jxz.begin(), jxz.end(), temploc) != jxz.end()) // Check for duplicate entries
                master->print_warning("Removed duplicate entry y=%f for [cross][xz]=%f\n", ycross, *it);
            else // Add to cross-list
            {
                jxz.push_back(temploc);
                master->print_message("Added XZ cross at y=%f (j=%i) for [cross][xz]=%f\n", ycross, temploc, *it);
            }

            if (std::find(jxzh.begin(), jxzh.end(), temploch) != jxzh.end()) // Check for duplicate entries
                master->print_warning("Removed duplicate entry yh=%f for [cross][xz]=%f\n", ycrossh, *it);
            else // Add to cross-list
            {
                jxzh.push_back(temploch);
                master->print_message("Added XZ cross at yh=%f (j=%i) for [cross][xz]=%f\n", ycrossh, temploch, *it);
            }
        }
    }
    
    // Find nearest full and half grid locations of yz cross-sections.
    for (std::vector<double>::iterator it=yz.begin(); it<yz.end(); ++it)
    {
        temploc  = (int) floor(*it/(grid->dx));
        temploch = (int) floor((*it+(grid->dx/2.))/(grid->dx));

        if (*it < 0 || *it > grid->xsize) // Check if cross location is inside domain
        {
            master->print_error("%f in [cross][yz] is outside domain\n", *it);
            ++nerror;
        }
        else
        {
            if (*it == grid->xsize) // Exception for full level when requesting domain size
                --temploc;

            // Find the corresponding index, make sure to handle MPI properly.
            double xcross = -1.;
            if (temploc / grid->imax == master->mpicoordx)
                xcross = grid->x[temploc % grid->imax + grid->igc];
            master->max(&xcross, 1);

            double xcrossh = -1.;
            if (temploch / grid->imax == master->mpicoordx)
                xcrossh = grid->x[temploch % grid->imax + grid->igc];
            master->max(&xcrossh, 1);

            if (std::find(ixz.begin(), ixz.end(), temploc) != ixz.end()) // Check for duplicate entries
                master->print_warning("Removed duplicate entry x=%f for [cross][yz]=%f\n", xcross, *it);
            else // Add to cross-list
            {
                ixz.push_back(temploc);
                master->print_message("Added YZ cross at x=%f (i=%i) for [cross][yz]=%f\n", xcross, temploc, *it);
            } 

            if (std::find(ixzh.begin(), ixzh.end(), temploch) != ixzh.end()) // Check for duplicate entries
                master->print_warning("Removed duplicate entry xh=%f for [cross][yz]=%f\n", xcrossh, *it);
            else // Add to cross-list
            {
                ixzh.push_back(temploch);
                master->print_message("Added YZ cross at xh=%f (i=%i) for [cross][yz]=%f\n", xcrossh, temploch, *it);
            } 
        }
    }

    // Find nearest full and half grid locations of xy cross-sections.
    for (std::vector<double>::iterator it=xy.begin(); it<xy.end(); ++it)
    {
        hoffset = 0;
        if (*it < 0 || *it > grid->zsize) // Check if cross location is inside domain
        {
            master->print_error("%f in [cross][xy] is outside domain\n", *it);
            ++nerror;
        }
        else
        {
            if (*it == grid->zsize) // Exception for domain top: use half level at domain top, full level below
            {
                temploc = grid->kmax-1;
                hoffset = 1;
            }
            else
            {
                for (int k=grid->kstart; k<grid->kend; k++) // Loop over height to find the nearest full level
                {
                    if ((*it >= grid->zh[k]) && (*it < grid->zh[k+1]))
                    {
                        temploc = k - grid->kgc;
                        if (*it >= grid->z[k]) // Add offset for half level
                            hoffset = 1;
                        break;
                    }
                }
            }

            if (std::find(kxy.begin(), kxy.end(), temploc) != kxy.end()) // Check for duplicate entries
                master->print_warning("Removed duplicate entry z=%f for [cross][xy]=%f\n", grid->z[temploc+grid->kgc],*it);
            else // Add to cross-list
            {
                kxy.push_back(temploc);
                master->print_message("Added XY cross at z=%f (k=%i) for [cross][xy]=%f\n", grid->z[temploc+grid->kgc],temploc,*it);
            } 

            if (std::find(kxyh.begin(), kxyh.end(), temploc+hoffset) != kxyh.end()) // Check for duplicate entries
                master->print_warning("Removed duplicate entry zh=%f for [cross][xy]=%f\n", grid->zh[temploc+hoffset+grid->kgc],*it);
            else // Add to cross-list
            {
                kxyh.push_back(temploc+hoffset);
                master->print_message("Added XY cross at zh=%f (k=%i) for [cross][xy]=%f\n", grid->zh[temploc+hoffset+grid->kgc],temploc+hoffset,*it);
            }  
        }
    }

    /* All classes (fields, thermo, boundary) have removed their cross-variables from
       crosslist by now. If it isnt empty, print warnings for invalid variables */
    if (crosslist.size() > 0)
    {
        for (std::vector<std::string>::const_iterator it=crosslist.begin(); it!=crosslist.end(); ++it)
            master->print_warning("field %s in [cross][crosslist] is illegal\n", it->c_str());
    } 

    if (nerror)
        throw 1;
}

unsigned long Cross::get_time_limit(unsigned long itime)
{
    if (swcross == "0")
        return Constants::ulhuge;

    unsigned long idtlim = isampletime - itime % isampletime;

    return idtlim;
}

std::string Cross::get_switch()
{
    return swcross;
}

std::vector<std::string>* Cross::get_crosslist()
{
    return &crosslist;
}

bool Cross::do_cross()
{
    if (swcross == "0")
        return false;

    if (model->timeloop->get_itime() % isampletime == 0)
        return true;
    else
        return false;
}

int Cross::cross_simple(double* restrict data, double* restrict tmp, std::string name)
{
    int nerror = 0;
    char filename[256];

    // loop over the index arrays to save all xz cross sections
    if (name == "v")
    {
        for (std::vector<int>::iterator it=jxzh.begin(); it<jxzh.end(); ++it)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", *it, model->timeloop->get_iotime());
            nerror += check_save(grid->save_xz_slice(data, tmp, filename, *it), filename);    
        }
    }
    else
    {
        for (std::vector<int>::iterator it=jxz.begin(); it<jxz.end(); ++it)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", *it, model->timeloop->get_iotime());
            nerror += check_save(grid->save_xz_slice(data, tmp, filename, *it), filename);    
        }
    }
    
    // loop over the index arrays to save all yz cross sections
    if (name == "u")
    {
        for (std::vector<int>::iterator it=ixzh.begin(); it<ixzh.end(); ++it)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "yz", *it, model->timeloop->get_iotime());
            nerror += check_save(grid->save_yz_slice(data, tmp, filename, *it), filename);    
        }
    }
    else
    {
        for (std::vector<int>::iterator it=ixz.begin(); it<ixz.end(); ++it)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "yz", *it, model->timeloop->get_iotime());
            nerror += check_save(grid->save_yz_slice(data, tmp, filename, *it), filename);    
        }
    }

    if (name == "w")
    {
        // loop over the index arrays to save all xy cross sections
        for (std::vector<int>::iterator it=kxyh.begin(); it<kxyh.end(); ++it)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", *it, model->timeloop->get_iotime());
            nerror += check_save(grid->save_xy_slice(data, tmp, filename, *it), filename);
        }
    }
    else
    {
        for (std::vector<int>::iterator it=kxy.begin(); it<kxy.end(); ++it)
        {
            std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", *it, model->timeloop->get_iotime());
            nerror += check_save(grid->save_xy_slice(data, tmp, filename, *it), filename);
        }
    }

    return nerror;
}

int Cross::cross_plane(double* restrict data, double* restrict tmp, std::string name)
{
    int nerror = 0;
    char filename[256];

    std::sprintf(filename, "%s.%s.%07d", name.c_str(), "xy", model->timeloop->get_iotime());
    nerror += check_save(grid->save_xy_slice(data, tmp, filename), filename);

    return nerror;
} 

int Cross::cross_lngrad(double* restrict a, double* restrict lngrad, double* restrict tmp, double* restrict dzi4, std::string name)
{
    using namespace Finite_difference::O4;

    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int jj3 = 3*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;
    const int kk3 = 3*grid->ijcells;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    int nerror = 0;
    char filename[256];

    // calculate the log of the gradient
    // bottom
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            lngrad[ijk] = std::log( Constants::dtiny +
                std::pow( ( cg0*(ci0*a[ijk-ii3] + ci1*a[ijk-ii2] + ci2*a[ijk-ii1] + ci3*a[ijk    ])
                          + cg1*(ci0*a[ijk-ii2] + ci1*a[ijk-ii1] + ci2*a[ijk    ] + ci3*a[ijk+ii1])
                          + cg2*(ci0*a[ijk-ii1] + ci1*a[ijk    ] + ci2*a[ijk+ii1] + ci3*a[ijk+ii2])
                          + cg3*(ci0*a[ijk    ] + ci1*a[ijk+ii1] + ci2*a[ijk+ii2] + ci3*a[ijk+ii3]) ) * cgi*dxi, 2.)

              + std::pow( ( cg0*(ci0*a[ijk-jj3] + ci1*a[ijk-jj2] + ci2*a[ijk-jj1] + ci3*a[ijk    ])
                          + cg1*(ci0*a[ijk-jj2] + ci1*a[ijk-jj1] + ci2*a[ijk    ] + ci3*a[ijk+jj1])
                          + cg2*(ci0*a[ijk-jj1] + ci1*a[ijk    ] + ci2*a[ijk+jj1] + ci3*a[ijk+jj2])
                          + cg3*(ci0*a[ijk    ] + ci1*a[ijk+jj1] + ci2*a[ijk+jj2] + ci3*a[ijk+jj3]) ) * cgi*dyi, 2.)

              + std::pow( ( cg0*(bi0*a[ijk-kk2] + bi1*a[ijk-kk1] + bi2*a[ijk    ] + bi3*a[ijk+kk1])
                          + cg1*(ci0*a[ijk-kk2] + ci1*a[ijk-kk1] + ci2*a[ijk    ] + ci3*a[ijk+kk1])
                          + cg2*(ci0*a[ijk-kk1] + ci1*a[ijk    ] + ci2*a[ijk+kk1] + ci3*a[ijk+kk2])
                          + cg3*(ci0*a[ijk    ] + ci1*a[ijk+kk1] + ci2*a[ijk+kk2] + ci3*a[ijk+kk3]) ) * dzi4[kstart], 2.) );
        }

    // interior
    for (int k=grid->kstart+1; k<grid->kend-1; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                lngrad[ijk] = std::log( Constants::dtiny +
                    std::pow( ( cg0*(ci0*a[ijk-ii3] + ci1*a[ijk-ii2] + ci2*a[ijk-ii1] + ci3*a[ijk    ])
                              + cg1*(ci0*a[ijk-ii2] + ci1*a[ijk-ii1] + ci2*a[ijk    ] + ci3*a[ijk+ii1])
                              + cg2*(ci0*a[ijk-ii1] + ci1*a[ijk    ] + ci2*a[ijk+ii1] + ci3*a[ijk+ii2])
                              + cg3*(ci0*a[ijk    ] + ci1*a[ijk+ii1] + ci2*a[ijk+ii2] + ci3*a[ijk+ii3]) ) * cgi*dxi, 2.)

                  + std::pow( ( cg0*(ci0*a[ijk-jj3] + ci1*a[ijk-jj2] + ci2*a[ijk-jj1] + ci3*a[ijk    ])
                              + cg1*(ci0*a[ijk-jj2] + ci1*a[ijk-jj1] + ci2*a[ijk    ] + ci3*a[ijk+jj1])
                              + cg2*(ci0*a[ijk-jj1] + ci1*a[ijk    ] + ci2*a[ijk+jj1] + ci3*a[ijk+jj2])
                              + cg3*(ci0*a[ijk    ] + ci1*a[ijk+jj1] + ci2*a[ijk+jj2] + ci3*a[ijk+jj3]) ) * cgi*dyi, 2.)

                  + std::pow( ( cg0*(ci0*a[ijk-kk3] + ci1*a[ijk-kk2] + ci2*a[ijk-kk1] + ci3*a[ijk    ])
                              + cg1*(ci0*a[ijk-kk2] + ci1*a[ijk-kk1] + ci2*a[ijk    ] + ci3*a[ijk+kk1])
                              + cg2*(ci0*a[ijk-kk1] + ci1*a[ijk    ] + ci2*a[ijk+kk1] + ci3*a[ijk+kk2])
                              + cg3*(ci0*a[ijk    ] + ci1*a[ijk+kk1] + ci2*a[ijk+kk2] + ci3*a[ijk+kk3]) ) * dzi4[k], 2.) );
            }

    // top
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            lngrad[ijk] = std::log(Constants::dtiny +
                std::pow( ( cg0*(ci0*a[ijk-ii3] + ci1*a[ijk-ii2] + ci2*a[ijk-ii1] + ci3*a[ijk    ])
                          + cg1*(ci0*a[ijk-ii2] + ci1*a[ijk-ii1] + ci2*a[ijk    ] + ci3*a[ijk+ii1])
                          + cg2*(ci0*a[ijk-ii1] + ci1*a[ijk    ] + ci2*a[ijk+ii1] + ci3*a[ijk+ii2])
                          + cg3*(ci0*a[ijk    ] + ci1*a[ijk+ii1] + ci2*a[ijk+ii2] + ci3*a[ijk+ii3]) ) * cgi*dxi, 2.)

              + std::pow( ( cg0*(ci0*a[ijk-jj3] + ci1*a[ijk-jj2] + ci2*a[ijk-jj1] + ci3*a[ijk    ])
                          + cg1*(ci0*a[ijk-jj2] + ci1*a[ijk-jj1] + ci2*a[ijk    ] + ci3*a[ijk+jj1])
                          + cg2*(ci0*a[ijk-jj1] + ci1*a[ijk    ] + ci2*a[ijk+jj1] + ci3*a[ijk+jj2])
                          + cg3*(ci0*a[ijk    ] + ci1*a[ijk+jj1] + ci2*a[ijk+jj2] + ci3*a[ijk+jj3]) ) * cgi*dyi, 2.)

              + std::pow( ( cg0*(ci0*a[ijk-kk3] + ci1*a[ijk-kk2] + ci2*a[ijk-kk1] + ci3*a[ijk    ])
                          + cg1*(ci0*a[ijk-kk2] + ci1*a[ijk-kk1] + ci2*a[ijk    ] + ci3*a[ijk+kk1])
                          + cg2*(ci0*a[ijk-kk1] + ci1*a[ijk    ] + ci2*a[ijk+kk1] + ci3*a[ijk+kk2])
                          + cg3*(ti0*a[ijk-kk1] + ti1*a[ijk    ] + ti2*a[ijk+kk1] + ti3*a[ijk+kk2]) ) * dzi4[kend-1], 2.) );
        }


    // loop over the index arrays to save all xz cross sections
    for (std::vector<int>::iterator it=jxz.begin(); it<jxz.end(); ++it)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", *it, model->timeloop->get_iotime());
        nerror += check_save(grid->save_xz_slice(lngrad, tmp, filename, *it),filename);
    }
    
    // loop over the index arrays to save all yz cross sections
    for (std::vector<int>::iterator it=ixz.begin(); it<ixz.end(); ++it)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "yz", *it, model->timeloop->get_iotime());
        nerror += check_save(grid->save_yz_slice(lngrad, tmp, filename, *it),filename);
    }

    // loop over the index arrays to save all xy cross sections
    for (std::vector<int>::iterator it=kxy.begin(); it<kxy.end(); ++it)
    {
        std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", *it, model->timeloop->get_iotime());
        nerror += check_save(grid->save_xy_slice(lngrad, tmp, filename, *it),filename);
    }

    return nerror;
}

int Cross::cross_path(double* restrict data, double* restrict tmp, double* restrict tmp1, std::string name)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    int nerror = 0;

    // Path is integrated in first full level, set to zero first
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj + kstart*kk;
            tmp[ijk] = 0.;
        }

    // Integrate with height
    for (int k=kstart; k<grid->kend; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk1 = i + j*jj + kstart*kk;
                const int ijk  = i + j*jj + k*kk;
                tmp[ijk1] += fields->rhoref[k] * data[ijk] * grid->dz[k];       
            }

    nerror += cross_plane(&tmp[kstart*kk], tmp1, name);

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
 * @param Direction Switch for bottom-up (Bottom_to_top) or top-down (Top_to_bottom) 
 * @param name String containing the output name of the cross-section 
 */
int Cross::cross_height_threshold(double* restrict data, double* restrict height, double* restrict tmp1, double* restrict z, double threshold, Direction direction, std::string name)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    int nerror = 0;

    // Set height to NetCDF fill value
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij = i + j*jj;
            height[ij] = NC_FILL_DOUBLE;
        }

    if(direction == Bottom_to_top) // Find lowest grid level where data > threshold
    {
        
        for (int j=grid->jstart; j<grid->jend; j++)
            for (int i=grid->istart; i<grid->iend; i++)
                for (int k=kstart; k<grid->kend; k++)
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
    else if(direction == Top_to_bottom) // Find highest grid level where data > threshold
    {
        for (int j=grid->jstart; j<grid->jend; j++)
            for (int i=grid->istart; i<grid->iend; i++)
                for (int k=grid->kend-1; k>grid->kstart-1; k--)
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

    nerror += cross_plane(height, tmp1, name);

    return nerror;
}
