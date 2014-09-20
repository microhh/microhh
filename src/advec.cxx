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
#include "advec.h"
#include "defines.h"
#include "constants.h"
#include "master.h"
#include "model.h"

#include "advec_2.h"
#include "advec_2int4.h"
#include "advec_4.h"
#include "advec_4m.h"

cadvec::cadvec(cmodel *modelin, cinput *inputin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  int nerror = 0;
  nerror += inputin->getItem(&cflmax, "advec", "cflmax", "", 1.);

  if(nerror)
    throw 1;
}

cadvec::~cadvec()
{
}

unsigned long cadvec::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim = (unsigned long) constants::dbig;

  return idtlim;
}

double cadvec::getcfl(double dt)
{
  double cfl;

  cfl = constants::dsmall;

  return cfl;
}

void cadvec::exec()
{
}

cadvec* cadvec::factory(cmaster *masterin, cinput *inputin, cmodel *modelin, const std::string swspatialorder)
{
  std::string swadvec;
  if(inputin->getItem(&swadvec, "advec", "swadvec", "", swspatialorder))
    throw 1;

  if(swadvec == "0")
    return new cadvec(modelin, inputin);
  else if(swadvec == "2")
    return new cadvec_2(modelin, inputin);
  else if(swadvec == "2int4")
    return new cadvec_2int4(modelin, inputin);
  else if(swadvec == "4")
    return new cadvec_4(modelin, inputin);
  else if(swadvec == "4m")
    return new cadvec_4m(modelin, inputin);
  else
  {
    masterin->printError("\"%s\" is an illegal value for swadvec\n", swadvec.c_str());
    throw 1;
  }
}
