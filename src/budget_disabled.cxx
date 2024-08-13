/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

template<typename TF>
Budget_disabled<TF>::Budget_disabled(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Input& inputin) :
    Budget<TF>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin)
{}

template<typename TF>
Budget_disabled<TF>::~Budget_disabled()
{}

template<typename TF>
void Budget_disabled<TF>::init()
{}

template<typename TF>
void Budget_disabled<TF>::create(Stats<TF>&)
{}

template<typename TF>
void Budget_disabled<TF>::exec_stats(Stats<TF>&)
{}


#ifdef FLOAT_SINGLE
template class Budget_disabled<float>;
#else
template class Budget_disabled<double>;
#endif
