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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "boundary_surface_lsm.h"
#include "boundary.h"

#include "master.h"
#include "input.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "soil_field3d.h"
#include "diff.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "master.h"
#include "cross.h"
#include "column.h"
#include "monin_obukhov.h"
#include "fast_math.h"
#include "netcdf_interface.h"
#include "radiation.h"
#include "microphys.h"

#include "boundary_surface_kernels.h"
#include "land_surface_kernels.h"
#include "soil_kernels.h"


#ifdef USECUDA
template<typename TF>
void Boundary_surface_lsm<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
    throw std::runtime_error("ERROR: boundary_surface_lsm is not (yet) implemented on the GPU!\n");
}

template<typename TF>
void Boundary_surface_lsm<TF>::exec_column(Column<TF>& column)
{
    throw std::runtime_error("ERROR: boundary_surface_lsm is not (yet) implemented on the GPU!\n");
}
#endif

template class Boundary_surface_lsm<double>;
template class Boundary_surface_lsm<float>;
