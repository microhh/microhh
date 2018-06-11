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

#include <iostream>
#include <cstdio>
#include "master.h"
#include "grid.h"
#include "fields.h"

#include "constants.h"
#include "microphys.h"
#include "microphys_2mom_warm.h"

template<typename TF>
Microphys_2mom_warm<TF>::Microphys_2mom_warm(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    swmicrophys = Microphys_type::Warm_2mom;
}

template<typename TF>
Microphys_2mom_warm<TF>::~Microphys_2mom_warm()
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::init()
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::create(Input& inputin, Data_block& data_block, Stats<TF>& stats, Cross<TF>& cross, Dump<TF>& dump)
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec()
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_stats(Stats<TF>& stats, std::string mask_name,
                                         Field3d<TF>& mask_field, Field3d<TF>& mask_fieldh, const double dt)
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
}

template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
bool Microphys_2mom_warm<TF>::has_mask(std::string name)
{
    return false;
}

template<typename TF>
void Microphys_2mom_warm<TF>::get_mask(Field3d<TF>& mfield, Field3d<TF>& mfieldh, Stats<TF>& stats, std::string mask_name)
{
}

template class Microphys_2mom_warm<double>;
template class Microphys_2mom_warm<float>;
