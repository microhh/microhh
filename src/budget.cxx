/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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


template<typename TF>
Budget<TF>::Budget(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin),
    thermo(thermoin), diff(diffin), advec(advecin), force(forcein)
{}

template<typename TF>
Budget<TF>::~Budget()
{}

template<typename TF>
std::shared_ptr<Budget<TF>> Budget<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Stats<TF>& statsin, Input& inputin)
{
    std::string swbudget = inputin.get_item<std::string>("budget", "swbudget", "", "0");

    // If the stats is disabled, also disable the budget stats.
    if (statsin.get_switch() == false)
        swbudget = "0";

    if (swbudget == "0")
        return std::make_shared<Budget_disabled<TF>>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin);
    else if (swbudget == "2")
        return std::make_shared<Budget_2<TF>>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin);
    else if (swbudget == "4")
        return std::make_shared<Budget_4<TF>>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin);
    else
    {
        std::string error_message = swbudget + " is an illegal value for swbudget";
        throw std::runtime_error(error_message);
    }
}


#ifdef FLOAT_SINGLE
template class Budget<float>;
#else
template class Budget<double>;
#endif
