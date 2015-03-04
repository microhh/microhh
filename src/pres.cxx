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

Pres::Pres(Model *modelin, Input *input)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  #ifdef USECUDA
  iplanf = 0;
  jplanf = 0;
  iplanb = 0;
  jplanb = 0;
  #endif
}

Pres::~Pres()
{
  #ifdef USECUDA
  cufftDestroy(iplanf);
  cufftDestroy(jplanf);
  cufftDestroy(iplanb);
  cufftDestroy(jplanb);
  #endif
}

void Pres::init()
{
}

void Pres::set_values()
{
}

void Pres::exec(double dt)
{
}

double Pres::checkDivergence()
{
  double divmax = 0.;
  return divmax;
}

Pres* Pres::factory(Master *masterin, Input *inputin, Model *modelin, const std::string swspatialorder)
{
  std::string swpres;
  if (inputin->getItem(&swpres, "pres", "swpres", "", swspatialorder))
    return 0;

  if (swpres == "0")
    return new Pres(modelin, inputin);
  else if (swpres == "2")
    return new Pres_2(modelin, inputin);
  else if (swpres == "4")
    return new Pres_4(modelin, inputin);
  else
  {
    masterin->print_error("\"%s\" is an illegal value for swpres\n", swpres.c_str());
    throw 1;
  }
}

void Pres::prepare_device()
{
}
