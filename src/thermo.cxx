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

void Thermo::exec()
{
}

void Thermo::exec_stats(Mask *f)
{
}

void Thermo::exec_cross()
{
}

void Thermo::exec_dump()
{
}

bool Thermo::checkThermoField(std::string name)
{
  return true;  // always returns error 
}

void Thermo::getThermoField(Field3d *field, Field3d *tmp, std::string name)
{
}

void Thermo::getBuoyancySurf(Field3d *bfield)
{
}

void Thermo::getBuoyancyFluxbot(Field3d *bfield)
{
}

std::string Thermo::getSwitch()
{
  return swthermo;
}

void Thermo::getProgVars(std::vector<std::string> *list)
{
}

void Thermo::get_mask(Field3d *mfield, Field3d *mfieldh, Mask *f)
{
}

Thermo* Thermo::factory(Master *masterin, Input *inputin, Model *modelin)
{
  std::string swthermo;
  if (inputin->getItem(&swthermo, "thermo", "swthermo", "", "0"))
    return 0;

  if (swthermo== "moist")
    return new Thermo_moist(modelin, inputin);
  else if (swthermo == "buoy")
    return new Thermo_buoy(modelin, inputin);
  else if (swthermo == "dry")
    return new Thermo_dry(modelin, inputin);
  else if (swthermo == "buoy_slope")
    return new Thermo_buoy_slope(modelin, inputin);
  else if (swthermo == "0")
    return new Thermo(modelin, inputin);
  else
  {
    masterin->print_error("\"%s\" is an illegal value for swthermo\n", swthermo.c_str());
    return 0;
  }
}

void Thermo::prepare_device()
{
}

void Thermo::clear_device()
{
}
