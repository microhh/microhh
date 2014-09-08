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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "model.h"

// thermo schemes
#include "thermo.h"
#include "thermo_buoy.h"
#include "thermo_buoy_slope.h"
#include "thermo_dry.h"
#include "thermo_moist.h"

cthermo::cthermo(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  swthermo = "0";
}

cthermo::~cthermo()
{
}

int cthermo::readinifile(cinput *inputin)
{
  return 0;
}

int cthermo::init()
{
  return 0;
}

int cthermo::create()
{
  return 0;
}

int cthermo::exec()
{
  return 0;
}

int cthermo::execstats(mask *f)
{
  return 0;
}

int cthermo::execcross()
{
  return 0;
}

int cthermo::checkthermofield(std::string name)
{
  return 1;  // always returns error 
}

int cthermo::getthermofield(cfield3d *field, cfield3d *tmp, std::string name)
{
  return 0;
}

int cthermo::getbuoyancysurf(cfield3d *bfield)
{
  return 0;
}

int cthermo::getbuoyancyfluxbot(cfield3d *bfield)
{
  return 0;
}

std::string cthermo::getsw()
{
  return swthermo;
}

int cthermo::getprogvars(std::vector<std::string> *list)
{
  return 0;
}

int cthermo::getmask(cfield3d *mfield, cfield3d *mfieldh, mask *f)
{
  return 0;
}

cthermo* cthermo::factory(cmaster *masterin, cinput *inputin, cmodel *modelin)
{
  std::string swthermo;
  if(inputin->getItem(&swthermo, "thermo", "swthermo", "", "0"))
    return 0;

  if(swthermo== "moist")
    return new cthermo_moist(modelin);
  else if(swthermo == "buoy")
    return new cthermo_buoy(modelin);
  else if(swthermo == "dry")
    return new cthermo_dry(modelin);
  else if(swthermo == "buoy_slope")
    return new cthermo_buoy_slope(modelin);
  else if(swthermo == "0")
    return new cthermo(modelin);
  else
  {
    masterin->printError("\"%s\" is an illegal value for swthermo\n", swthermo.c_str());
    return 0;
  }
}
