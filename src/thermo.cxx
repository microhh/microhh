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
#include "thermo.h"
#include "defines.h"
#include "model.h"

cthermo::cthermo(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  swthermo = "off";
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

int cthermo::execstats()
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

std::string cthermo::getname()
{
  return swthermo;
}

int cthermo::getprogvars(std::vector<std::string> *list)
{
  return 0;
}
