/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include "radiation.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "radiation_rrtmgp.h"

template<typename TF>
Radiation_rrtmgp<TF>::Radiation_rrtmgp(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
	Radiation<TF>(masterin, gridin, fieldsin, inputin)
{
	swradiation = "rrtmgp";
}

template<typename TF>
Radiation_rrtmgp<TF>::~Radiation_rrtmgp()
{
}

template<typename TF>
void Radiation_rrtmgp<TF>::init()
{
}

template<typename TF>
void Radiation_rrtmgp<TF>::create(Thermo<TF>& thermo, Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
}

template<typename TF>
void Radiation_rrtmgp<TF>::exec(Thermo<TF>& thermo, const double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
}

template class Radiation_rrtmgp<double>;
template class Radiation_rrtmgp<float>;
