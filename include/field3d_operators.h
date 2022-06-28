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

#ifndef FIELD3D_OPERATORS_H
#define FIELD3D_OPERATORS_H

#include "field3d.h"

class Master;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Field3d;

template<typename TF>
class Field3d_operators
{
    public:
        Field3d_operators(Master&, Grid<TF>&, Fields<TF>&);
        ~Field3d_operators();

        void calc_mean_profile(TF* const, const TF* const); // Calculate mean profile into fld_mean
        void calc_mean_profile_nogc(TF* const, const TF* const, bool); // Calculate mean profile into fld_mean
        void subtract_mean_profile(TF* const, const TF* const); // Calculate mean profile into fld_mean
        TF calc_mean_2d(const TF* const); // Calculate mean from 2D field
        TF calc_mean(const TF* const); // Calculate volume-weighted mean.
        TF calc_max(const TF* const);

        #ifdef USECUDA
        void calc_mean_profile_g(TF* const, const TF* const); // Calculate mean profile into fld_mean
        TF calc_mean_2d_g(const TF* const); // Calculate mean from 2D (xy) field.
        TF calc_mean_g(const TF* const); // Calculate volume-weighted mean.
        TF calc_max_g(const TF* const);
        #endif

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
};
#endif
