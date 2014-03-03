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
#include "stats.h"
#include "defines.h"
#include "model.h"

cstats::cstats(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;
}

cstats::~cstats()
{
}

int cstats::readinifile(cinput *inputin)
{
  return 0;
}

int cstats::init(double ifactor)
{
  return 0;
}

int cstats::create(int n)
{
  return 0;
}

unsigned long cstats::gettimelim(unsigned long itime)
{
  return ulhuge;
}

int cstats::exec(int iteration, double time, unsigned long itime)
{
  return 0;
}

