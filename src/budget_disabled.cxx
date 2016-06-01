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

#include "budget_disabled.h"

Budget_disabled::Budget_disabled(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Diff* diffin, Advec* advecin, Force* forcein, Stats* statsin) :
    Budget(inputin, masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, statsin)
{
}

Budget_disabled::~Budget_disabled()
{
}

void Budget_disabled::init()
{
}

void Budget_disabled::create()
{
}

void Budget_disabled::exec_stats(Mask* m)
{
}
