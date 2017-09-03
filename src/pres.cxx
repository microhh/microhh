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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "model.h"

#include "pres.h"
#include "pres_2.h"
//#include "pres_4.h"

template<typename TF>
Pres<TF>::Pres(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
#ifdef USECUDA
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
double Pres<TF>::check_divergence()
{
    double divmax = 0.;
    return divmax;
}

template<typename TF>
std::shared_ptr<Pres<TF>> Pres<TF>::factory(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Input& inputin, const std::string swspatialorder)
{
    std::string swpres = inputin.get_item<std::string>("pres", "swpres", "", swspatialorder);

    if (swpres == "0")
        return std::make_shared<Pres<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swpres == "2")
        return std::make_shared<Pres_2<TF>>(masterin, gridin, fieldsin, inputin);
    // else if (swpres == "4")
    //     return new Pres_4(modelin, inputin);
    else
    {
        masterin.print_error("\"%s\" is an illegal value for swpres\n", swpres.c_str());
        throw 1;
    }
}

template<typename TF>
void Pres<TF>::prepare_device()
{
}

template class Pres<double>;
template class Pres<float>;
