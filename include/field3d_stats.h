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

#ifndef FIELD3D_STATS
#define FIELD3D_STATS

#include "field3d.h"

class Master;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Field3d;

template<typename TF>
class Field3d_stats
{
    public:
        // Functions
        Field3d_stats(Master&, Grid<TF>&, Fields<TF>&);
        ~Field3d_stats();


        void calc_mean_profile(Field3d<TF>*);   ///< Calculate mean profile into fld_mean
        TF calc_mean(Field3d<TF>*);             ///< Calculate volume weighted total mean
        TF calc_max(Field3d<TF>*);              ///< Calculate maximum value


    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
};
#endif

