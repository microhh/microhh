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

#include "input.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"

#include "budget.h"
#include "budget_disabled.h"

Budget::Budget(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Stats* statsin):
    master(*masterin),
    grid  (*gridin  ),
    fields(*fieldsin),
    thermo(*thermoin),
    stats (*statsin )
{
}

Budget::~Budget()
{
}

Budget* Budget::factory(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Stats* statsin)
{
    std::string swbudget;
    if (inputin->get_item(&swbudget, "budget", "swbudget", "", "0"))
        throw 1;

    if (swbudget == "0")
        return new Budget_disabled(inputin, masterin, gridin, fieldsin, thermoin, statsin);
    else
    {
        masterin->print_error("\"%s\" is an illegal value for swbudget\n", swbudget.c_str());
        throw 1;
    }
}
