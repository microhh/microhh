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
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "defines.h"

cthermo::cthermo(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  swthermo = "off";
}

cthermo::~cthermo()
{
}

int cthermo::readinifile(cinput *inputin)
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

int cthermo::getbuoyancysurf(cfield3d *bfield)
{
  return 0;
}

int cthermo::getbuoyancyfluxbot(cfield3d *bfield)
{
  return 0;
}

int cthermo::getbuoyancy(cfield3d *bfield, cfield3d *tmp)
{
  return 0;
}

std::string cthermo::getname()
{
  return swthermo;
}
