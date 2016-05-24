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
#include "defines.h"
#include "constants.h"
#include "master.h"
#include "model.h"

#include "advec.h"
#include "advec_disabled.h"
#include "advec_2.h"
#include "advec_2i4.h"
#include "advec_4.h"
#include "advec_4m.h"

Advec::Advec(Model* modelin, Input* inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    int nerror = 0;
    nerror += inputin->get_item(&cflmax, "advec", "cflmax", "", 1.);

    swadvec = "0";

    if (nerror)
        throw 1;
}

Advec::~Advec()
{
}

Advec* Advec::factory(Master* masterin, Input* inputin, Model* modelin, const std::string swspatialorder)
{
    std::string swadvec;
    if (inputin->get_item(&swadvec, "advec", "swadvec", "", swspatialorder))
        throw 1;

    if (swadvec == "0")
        return new Advec_disabled(modelin, inputin);
    else if (swadvec == "2")
        return new Advec_2(modelin, inputin);
    else if (swadvec == "2i4")
        return new Advec_2i4(modelin, inputin);
    else if (swadvec == "4")
        return new Advec_4(modelin, inputin);
    else if (swadvec == "4m")
        return new Advec_4m(modelin, inputin);
    else
    {
        masterin->print_error("\"%s\" is an illegal value for swadvec\n", swadvec.c_str());
        throw 1;
    }
}

std::string Advec::get_switch()
{
    return swadvec;
}

const double Advec::cflmin = 1.E-5;
