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
#include "master.h"
#include "diff.h"
#include "defines.h"
#include "model.h"

cdiff::cdiff(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  mpi    = model->mpi;
}

cdiff::~cdiff()
{
}

int cdiff::readinifile(cinput *inputin)
{
  int nerror = 0;

  nerror += inputin->getItem(&dnmax, "diff", "dnmax", "", 0.4);

  return nerror;
}

unsigned long cdiff::gettimelim(unsigned long idtlim, double dt)
{
  idtlim = (unsigned long) dbig;

  return idtlim;
}

int cdiff::setvalues()
{
  return 0;
}

double cdiff::getdn(double dt)
{
  double dn;

  dn = dsmall;

  return dn;
}

int cdiff::execvisc()
{
  return 0;
}

int cdiff::exec()
{
  return 0;
}

