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
#include <iostream>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "field3d.h"
#include "defines.h"

cfield3d::cfield3d(cgrid *gridin, cmaster *masterin, std::string namein, std::string longnamein, std::string unitin)
{
  grid     = gridin;
  name     = namein;
  longname = longnamein;
  unit     = unitin;
  master   = masterin;

  // initialize the pointer at 0
  data = 0;
  databot = 0;
  datatop = 0;
  datamean = 0;
  datagradbot = 0;
  datagradtop = 0;
  datafluxbot = 0;
  datafluxtop = 0;
}

cfield3d::~cfield3d()
{
  delete[] data;
  delete[] databot;
  delete[] datatop;
  delete[] datamean;
  delete[] datagradbot;
  delete[] datagradtop;
  delete[] datafluxbot;
  delete[] datafluxtop;
}

int cfield3d::init()
{
  // allocate the memory
  master->printMessage("Allocating %d bytes of memory for %s\n", grid->ncells*(int)sizeof(double), name.c_str());
  data = new double[grid->ncells];

  // allocate the boundary cells
  databot = new double[grid->icells*grid->jcells];
  datatop = new double[grid->icells*grid->jcells];
  datamean = new double[grid->kcells];
  datagradbot = new double[grid->icells*grid->jcells];
  datagradtop = new double[grid->icells*grid->jcells];
  datafluxbot = new double[grid->icells*grid->jcells];
  datafluxtop = new double[grid->icells*grid->jcells];

  // set all values to zero
  for(int n=0; n<grid->ncells; n++)
    data[n] = 0.;

  for(int n=0; n<grid->kcells; n++)
    datamean[n] = 0.;

  for(int n=0; n<grid->icells*grid->jcells; n++)
  {
    databot    [n] = 0.;
    datatop    [n] = 0.;
    datagradbot[n] = 0.;
    datagradtop[n] = 0.;
    datafluxbot[n] = 0.;
    datafluxtop[n] = 0.;
  }

  return 0;
}

/*
int cfield3d::checkfornan()
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;
  int nerror=0;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double cfl = 0.;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        if (data[ijk]!=data[ijk])
        {
          nerror++;
          std::printf("Cell %d %d %d of %s is %f\n", i,j,k, name.c_str(),data[ijk]);
        }
      }

  return nerror;
}
*/
