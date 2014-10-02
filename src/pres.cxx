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
#include <fftw3.h>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "model.h"

#include "pres.h"
#include "pres_2.h"
#include "pres_4.h"

cpres::cpres(cmodel *modelin, cinput *input)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;
}

cpres::~cpres()
{
}

void cpres::init()
{
}

void cpres::setvalues()
{
}

void cpres::exec(double dt)
{
}

double cpres::check()
{
  double divmax = 0.;
  return divmax;
}

cpres* cpres::factory(cmaster *masterin, cinput *inputin, cmodel *modelin, const std::string swspatialorder)
{
  std::string swpres;
  if(inputin->getItem(&swpres, "pres", "swpres", "", swspatialorder))
    return 0;

  if(swpres == "0")
    return new cpres(modelin, inputin);
  else if(swpres == "2")
    return new cpres_2(modelin, inputin);
  else if(swpres == "4")
    return new cpres_4(modelin, inputin);
  else
  {
    masterin->printError("\"%s\" is an illegal value for swpres\n", swpres.c_str());
    throw 1;
  }
}

int cpres::prepareDevice()
{
  return 0;
}
