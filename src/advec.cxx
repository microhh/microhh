/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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
#include "advec.h"
#include "defines.h"
#include "model.h"

#include "advec_2.h"
#include "advec_2int4.h"
#include "advec_4.h"
#include "advec_4m.h"

cadvec::cadvec(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;
}

cadvec::~cadvec()
{
}

int cadvec::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&cflmax, "advec", "cflmax", "", 1.);

  return nerror;
}

unsigned long cadvec::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim = (unsigned long) dbig;

  return idtlim;
}

double cadvec::getcfl(double dt)
{
  double cfl;

  cfl = dsmall;

  return cfl;
}

int cadvec::exec()
{
  return 0;
}

cadvec* cadvec::factory(cmodel *modelin, std::string swadvec)
{
  // check the advection scheme
  if(swadvec == "0")
    return new cadvec(modelin);
  else if(swadvec == "2")
    return new cadvec_2(modelin);
  else if(swadvec == "2int4")
    return new cadvec_2int4(modelin);
  else if(swadvec == "4")
    return new cadvec_4(modelin);
  else if(swadvec == "4m")
    return new cadvec_4m(modelin);
  else
    return 0;
}
