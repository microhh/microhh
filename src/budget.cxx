/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "diff.h"
#include "stats.h"

#include "budget.h"
#include "budget_disabled.h"
#include "budget_2.h"
#include "budget_4.h"

Budget::Budget(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Diff* diffin, Advec* advecin, Force* forcein, Stats* statsin):
    master(*masterin),
    grid  (*gridin  ),
    fields(*fieldsin),
    thermo(*thermoin),
    diff  (*diffin  ),
    advec (*advecin ),
    force (*forcein ),
    stats (*statsin )
{
}

Budget::~Budget()
{
}

Budget* Budget::factory(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Diff* diffin, Advec* advecin, Force* forcein, Stats* statsin)
{
    std::string swbudget;
    if (inputin->get_item(&swbudget, "budget", "swbudget", "", "0"))
        throw 1;

    // If the stats is disabled, also disable the budget stats.
    if (statsin->get_switch() == "0")
        swbudget = "0";

    if (swbudget == "0")
        return new Budget_disabled(inputin, masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, statsin);
    else if (swbudget == "2")
        return new Budget_2(inputin, masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, statsin);
    else if (swbudget == "4")
        return new Budget_4(inputin, masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, statsin);
    else
    {
        masterin->print_error("\"%s\" is an illegal value for swbudget\n", swbudget.c_str());
        throw 1;
    }
}
