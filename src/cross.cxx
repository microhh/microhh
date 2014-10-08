/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
#include "fd.h"
#include "model.h"
#include "thermo.h"
#include "timeloop.h"
#include <netcdfcpp.h>

Cross::Cross(Model *modelin, Input *inputin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  // optional, by default switch cross off
  int nerror = 0;
  nerror += inputin->getItem(&swcross, "cross", "swcross", "", "0");

  if(swcross == "1")
  {
    // get the time at which the cross sections are triggered
    nerror += inputin->getItem(&sampletime, "cross", "sampletime", "");

    // get the list of indices at which to take cross sections
    nerror += inputin->getList(&xz, "cross", "xz", "");
    nerror += inputin->getList(&xy, "cross", "xy", "");
  }

  if(nerror)
    throw 1;
}

Cross::~Cross()
{
}

// check whether saving the slice was successful and print appropriate message
int Cross::checkSave(int error, char * filename)
{
  master->printMessage("Saving \"%s\" ... ", filename);
  if(error == 0)
  {
    master->printMessage("OK\n");
    return 0;
  }
  else
  {
    master->printMessage("FAILED\n");
    return 1;
  }
}

void Cross::init(double ifactor)
{
  if(swcross == "0")
    return;

  isampletime = (unsigned long)(ifactor * sampletime);
}

void Cross::create()
{  
  int nerror = 0;
  int temploc, temploch, hoffset;

  // Find nearest full and half grid locations of xz cross-sections.
  for(std::vector<double>::iterator it=xz.begin(); it<xz.end(); ++it)
  {
    temploc  = (int) floor(*it/(grid->dy));
    temploch = (int) floor((*it+(grid->dy/2.))/(grid->dy));

    if(*it < 0 || *it > grid->ysize) // Check if cross location is inside domain
    {
      master->printError("ERROR %f in [cross][xz] is outside domain\n", *it);
      ++nerror;
    }
    else
    {
      if(*it == grid->ysize) // Exception for full level when requesting domain size
        --temploc;

      if(std::find(jxz.begin(), jxz.end(), temploc) != jxz.end()) // Check for duplicate entries
        master->printWarning("removed duplicate entry y=%f for [cross][xz]=%f\n", grid->y[temploc+grid->jgc],*it);
      else // Add to cross-list
      {
        jxz.push_back(temploc);
        master->printMessage("Addex XZ cross at y=%f (j=%i) for [cross][xz]=%f\n", grid->y[temploc+grid->jgc],temploc,*it);
      } 

      if(std::find(jxzh.begin(), jxzh.end(), temploch) != jxzh.end()) // Check for duplicate entries
        master->printWarning("removed duplicate entry yh=%f for [cross][xz]=%f\n", grid->yh[temploch+grid->jgc],*it);
      else // Add to cross-list
      {
        jxzh.push_back(temploch);
        master->printMessage("Addex XZ cross at yh=%f (j=%i) for [cross][xz]=%f\n", grid->yh[temploch+grid->jgc],temploch,*it);
      } 
    }
  }

  // Find nearest full and half grid locations of xy cross-sections.
  for(std::vector<double>::iterator it=xy.begin(); it<xy.end(); ++it)
  {
    hoffset = 0;
    if(*it < 0 || *it > grid->zsize) // Check if cross location is inside domain
    {
      master->printError("ERROR %f in [cross][xy] is outside domain\n", *it);
      ++nerror;
    }
    else
    {
      if(*it == grid->zsize) // Exception for domain top: use half level at domain top, full level below
      {
          temploc = grid->kmax-1;
          hoffset = 1;
      }
      else
      {
        for(int k=grid->kstart; k<grid->kend; k++) // Loop over height to find the nearest full level
        {
          if((*it >= grid->zh[k]) && (*it < grid->zh[k+1]))
          {
            temploc = k - grid->kgc;
            if(*it >= grid->z[k]) // Add offset for half level
              hoffset = 1;
            break;
          }
        }
      }

      if(std::find(kxy.begin(), kxy.end(), temploc) != kxy.end()) // Check for duplicate entries
        master->printWarning("removed duplicate entry z=%f for [cross][xy]=%f\n", grid->z[temploc+grid->kgc],*it);
      else // Add to cross-list
      {
        kxy.push_back(temploc);
        master->printMessage("Addex XY cross at z=%f (k=%i) for [cross][xy]=%f\n", grid->z[temploc+grid->kgc],temploc,*it);
      } 
 
      if(std::find(kxyh.begin(), kxyh.end(), temploc+hoffset) != kxyh.end()) // Check for duplicate entries
        master->printWarning("removed duplicate entry zh=%f for [cross][xy]=%f\n", grid->zh[temploc+hoffset+grid->kgc],*it);
      else // Add to cross-list
      {
        kxyh.push_back(temploc+hoffset);
        master->printMessage("Addex XY cross at zh=%f (k=%i) for [cross][xy]=%f\n", grid->zh[temploc+hoffset+grid->kgc],temploc+hoffset,*it);
      }  
    }
  }

  if(nerror)
    throw 1;
}

unsigned long Cross::getTimeLimit(unsigned long itime)
{
  if(swcross == "0")
    return constants::ulhuge;

  unsigned long idtlim = isampletime - itime % isampletime;

  return idtlim;
}

bool Cross::doCross()
{
  if(swcross == "0")
    return false;

  if(model->timeloop->get_itime() % isampletime == 0)
    return true;
  else
    return false;
}

int Cross::crossSimple(double * restrict data, double * restrict tmp, std::string name)
{
  int nerror = 0;
  char filename[256];

  // loop over the index arrays to save all xz cross sections
  if(name == "v")
  {
    for(std::vector<int>::iterator it=jxzh.begin(); it<jxzh.end(); ++it)
    {
      std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", *it, model->timeloop->get_iotime());
      nerror += checkSave(grid->savexzSlice(data, tmp, filename, *it), filename);    
    }
  }
  else
  {
    for(std::vector<int>::iterator it=jxz.begin(); it<jxz.end(); ++it)
    {
      std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", *it, model->timeloop->get_iotime());
      nerror += checkSave(grid->savexzSlice(data, tmp, filename, *it), filename);    
    }
  }

  if(name == "w")
  {
    // loop over the index arrays to save all xy cross sections
    for(std::vector<int>::iterator it=kxyh.begin(); it<kxyh.end(); ++it)
    {
      std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", *it, model->timeloop->get_iotime());
      nerror += checkSave(grid->savexySlice(data, tmp, filename, *it), filename);
    }
  }
  else
  {
    for(std::vector<int>::iterator it=kxy.begin(); it<kxy.end(); ++it)
    {
      std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", *it, model->timeloop->get_iotime());
      nerror += checkSave(grid->savexySlice(data, tmp, filename, *it), filename);
    }
  }

  return nerror;
}

int Cross::crossPlane(double * restrict data, double * restrict tmp, std::string name)
{
  int nerror = 0;
  char filename[256];

  std::sprintf(filename, "%s.%s.%07d", name.c_str(), "xy", model->timeloop->get_iotime());
  nerror += checkSave(grid->savexySlice(data, tmp, filename),filename);

  return nerror;
} 

int Cross::crossLngrad(double * restrict a, double * restrict lngrad, double * restrict tmp, double * restrict dzi4, std::string name)
{
  using namespace fd::o4;

  int ijk,ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  int kstart,kend;
  int nerror = 0;
  char filename[256];

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
  kk3 = 3*grid->ijcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  double dxi,dyi;
  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // calculate the log of the gradient
  // bottom
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk  = i + j*jj1 + kstart*kk1;
      lngrad[ijk] = std::log( constants::dtiny +
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
  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        lngrad[ijk] = std::log( constants::dtiny +
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
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk  = i + j*jj1 + (kend-1)*kk1;
      lngrad[ijk] = std::log(constants::dtiny +
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
  for(std::vector<int>::iterator it=jxz.begin(); it<jxz.end(); ++it)
  {
    std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xz", *it, model->timeloop->get_iotime());
    nerror += checkSave(grid->savexzSlice(lngrad, tmp, filename, *it),filename);
  }

  // loop over the index arrays to save all xy cross sections
  for(std::vector<int>::iterator it=kxy.begin(); it<kxy.end(); ++it)
  {
    std::sprintf(filename, "%s.%s.%05d.%07d", name.c_str(), "xy", *it, model->timeloop->get_iotime());
    nerror += checkSave(grid->savexySlice(lngrad, tmp, filename, *it),filename);
  }

  return nerror;
}

int Cross::crossPath(double * restrict data, double * restrict tmp, double * restrict tmp1, std::string name)
{
  int ijk,ijk1,jj,kk;
  jj = grid->icells;
  kk = grid->ijcells;
  int kstart = grid->kstart;
  int nerror = 0;
  char filename[256];

  // Path is integrated in first full level, set to zero first
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj + kstart*kk;
      tmp[ijk] = 0.;
    }

  // Integrate with height
  for(int k=kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk1 = i + j*jj + kstart*kk;
        ijk  = i + j*jj + k*kk;
        tmp[ijk1] += fields->rhoref[k] * data[ijk] * grid->dz[k];       
      }

  std::sprintf(filename, "%s.%s.%07d", name.c_str(), "xy", model->timeloop->get_iotime());
  nerror += checkSave(grid->savexySlice(&tmp[kstart*kk], tmp1, filename),filename);

  return nerror;
}
