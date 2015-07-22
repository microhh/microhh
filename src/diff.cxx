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
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "constants.h"
#include "model.h"

// diffusion schemes
#include "diff.h"
#include "diff_disabled.h"
#include "diff_2.h"
#include "diff_4.h"
#include "diff_smag2.h"

Diff::Diff(Model* modelin, Input* inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    swdiff = "0";

    int nerror = 0;
    nerror += inputin->get_item(&dnmax, "diff", "dnmax", "", 0.4);

    if (nerror)
        throw 1;
}

Diff::~Diff()
{
}

Diff* Diff::factory(Master* masterin, Input* inputin, Model* modelin, const std::string swspatialorder)
{
    std::string swdiff;
    std::string swboundary;

    int nerror = 0;
    nerror += inputin->get_item(&swdiff, "diff", "swdiff", "", swspatialorder);
    // load the boundary switch as well in order to be able to check whether the surface model is used
    nerror += inputin->get_item(&swboundary, "boundary", "swboundary", "", "default");
    if (nerror)
        return 0;

    if (swdiff == "0")
        return new Diff_disabled(modelin, inputin);
    else if (swdiff == "2")
        return new Diff_2(modelin, inputin);
    else if (swdiff == "4")
        return new Diff_4(modelin, inputin);
    else if (swdiff == "smag2")
    {
        // the subgrid model requires a surface model because of the MO matching at first level
        if ((swboundary == "surface") || (swboundary == "bulk"))
            return new Diff_smag_2(modelin, inputin);
        else
        {
            masterin->print_error("swdiff == \"smag2\" requires swboundary == \"surface\"\n");
            return 0;
        }
    }
    else
    {
        masterin->print_error("\"%s\" is an illegal value for swdiff\n", swdiff.c_str());
        return 0;
    }
}

std::string Diff::get_name()
{
    return swdiff;
}
