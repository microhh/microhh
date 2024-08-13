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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iostream>

#include <fftw3.h>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"

#include "pres.h"
#include "pres_disabled.h"
#include "pres_2.h"
#include "pres_4.h"

template<typename TF>
Pres<TF>::Pres(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, FFT<TF>& fftin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), fft(fftin),
    field3d_operators(master, grid, fields)
{
    #ifdef USECUDA
    force_FFT_per_slice = inputin.get_item<bool>("pres", "sw_fft_per_slice", "", false);

    iplanf = 0;
    jplanf = 0;
    iplanb = 0;
    jplanb = 0;
    #endif
}

template<typename TF>
Pres<TF>::~Pres()
{
    #ifdef USECUDA
    cufftDestroy(iplanf);
    cufftDestroy(jplanf);
    cufftDestroy(iplanb);
    cufftDestroy(jplanb);
    #endif
}

template<typename TF>
std::shared_ptr<Pres<TF>> Pres<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, FFT<TF>& fftin, Input& inputin)
{
    std::string swspatialorder = (gridin.get_spatial_order() == Grid_order::Second) ? "2" : "4";
    std::string swpres = inputin.get_item<std::string>("pres", "swpres", "", swspatialorder);

    if (swpres == "0")
        return std::make_shared<Pres_disabled<TF>>(masterin, gridin, fieldsin, fftin, inputin);
    else if (swpres == "2")
        return std::make_shared<Pres_2<TF>>(masterin, gridin, fieldsin, fftin, inputin);
    else if (swpres == "4")
        return std::make_shared<Pres_4<TF>>(masterin, gridin, fieldsin, fftin, inputin);
    else
    {
        std::string msg = swpres + " is an illegal value for swpres";
        throw std::runtime_error(msg);
    }
}


#ifdef FLOAT_SINGLE
template class Pres<float>;
#else
template class Pres<double>;
#endif
