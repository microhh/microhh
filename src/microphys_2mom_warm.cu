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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "boundary_cyclic.h"
#include "thermo_moist_functions.h"
#include "constants.h"

#include "microphys.h"
#include "microphys_2mom_warm.h"

using namespace Constants;
using namespace Thermo_moist_functions;

#ifdef USECUDA
template<typename TF>
void Microphys_2mom_warm<TF>::exec(Thermo<TF>& thermo, const double dt)
{
}
#endif

#ifdef USECUDA
template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}
#endif

template class Microphys_2mom_warm<double>;
template class Microphys_2mom_warm<float>;
