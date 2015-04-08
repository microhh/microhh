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

Diff_disabled::Diff_disabled(Model* modelin, Input* inputin) : Diff(modelin, inputin)
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

Diff_disabled::~Diff_disabled()
{
}

unsigned long Diff_disabled::get_time_limit(const unsigned long idtlim, const double dt)
{
    return Constants::ulhuge;
}

double Diff_disabled::get_dn(const double dt)
{
    return Constants::dsmall;
}
