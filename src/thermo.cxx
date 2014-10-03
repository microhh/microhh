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

Thermo::Thermo(Model *modelin, Input *inputin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  swthermo = "0";
}

Thermo::~Thermo()
{
}

void Thermo::init()
{
}

void Thermo::create(Input *inputin)
{
}

int Thermo::exec()
{
  return 0;
}

int Thermo::execstats(mask *f)
{
  return 0;
}

void Thermo::execcross()
{
}

int Thermo::checkthermofield(std::string name)
{
  return 1;  // always returns error 
}

int Thermo::getthermofield(Field3d *field, Field3d *tmp, std::string name)
{
  return 0;
}

int Thermo::getbuoyancysurf(Field3d *bfield)
{
  return 0;
}

int Thermo::getbuoyancyfluxbot(Field3d *bfield)
{
  return 0;
}

std::string Thermo::getsw()
{
  return swthermo;
}

int Thermo::getprogvars(std::vector<std::string> *list)
{
  return 0;
}

int Thermo::getmask(Field3d *mfield, Field3d *mfieldh, mask *f)
{
  return 0;
}

Thermo* Thermo::factory(Master *masterin, Input *inputin, Model *modelin)
{
  std::string swthermo;
  if(inputin->getItem(&swthermo, "thermo", "swthermo", "", "0"))
    return 0;

  if(swthermo== "moist")
    return new Thermo_moist(modelin, inputin);
  else if(swthermo == "buoy")
    return new Thermo_buoy(modelin, inputin);
  else if(swthermo == "dry")
    return new Thermo_dry(modelin, inputin);
  else if(swthermo == "buoy_slope")
    return new Thermo_buoy_slope(modelin, inputin);
  else if(swthermo == "0")
    return new Thermo(modelin, inputin);
  else
  {
    masterin->printError("\"%s\" is an illegal value for swthermo\n", swthermo.c_str());
    return 0;
  }
}

int Thermo::prepareDevice()
{
  return 0;
}

int Thermo::clearDevice()
{
  return 0;
}
