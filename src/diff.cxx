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
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "model.h"

// diffusion schemes
#include "diff.h"
#include "diff_2.h"
#include "diff_4.h"
#include "diff_les2s.h"

cdiff::cdiff(cmodel *modelin, cinput *inputin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  swdiff = "0";

  int nerror = 0;
  nerror += inputin->getItem(&dnmax, "diff", "dnmax", "", 0.4);
 
  if(nerror)
    throw 1;
}

cdiff::~cdiff()
{
}

unsigned long cdiff::gettimelim(unsigned long idtlim, double dt)
{
  idtlim = (unsigned long) dbig;

  return idtlim;
}

int cdiff::setvalues()
{
  return 0;
}

double cdiff::getdn(double dt)
{
  double dn;

  dn = dsmall;

  return dn;
}

int cdiff::execvisc()
{
  return 0;
}

int cdiff::exec()
{
  return 0;
}

std::string cdiff::getname()
{
  return swdiff;
}

cdiff* cdiff::factory(cmaster *masterin, cinput *inputin, cmodel *modelin, const std::string swspatialorder)
{
  std::string swdiff;
  std::string swboundary;

  int nerror = 0;
  nerror += inputin->getItem(&swdiff, "diff", "swdiff", "", swspatialorder);
  // load the boundary switch as well in order to be able to check whether the surface model is used
  nerror += inputin->getItem(&swboundary, "boundary", "swboundary", "", "default");
  if(nerror)
    return 0;

  if(swdiff == "0")
    return new cdiff(modelin, inputin);
  else if(swdiff == "2")
    return new cdiff_2(modelin, inputin);
  else if(swdiff == "4")
    return new cdiff_4(modelin, inputin);
  else if(swdiff == "les2s")
  {
    // the subgrid model requires a surface model because of the MO matching at first level
    if(swboundary != "surface")
    {
      masterin->printError("swdiff == \"les2s\" requires swboundary == \"surface\"\n");
      return 0;
    }
    return new cdiff_les2s(modelin, inputin);
  }
  else
  {
    masterin->printError("\"%s\" is an illegal value for swdiff\n", swdiff.c_str());
    return 0;
  }
}
