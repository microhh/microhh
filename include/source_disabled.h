/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#ifndef SOURCE_DISABLED_H
#define SOURCE_DISABLED_H

#include "source.h"

class Input;

template<typename TF>
class Source_disabled : public Source<TF>
{
    public:
        Source_disabled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Source_disabled();

        void init();
        void create(Input&, Netcdf_handle&);
        void exec();
        void update_time_dependent(Timeloop<TF>&);

        #ifdef USECUDA
        void prepare_device();
        #endif

    private:
        using Source<TF>::master;
        using Source<TF>::grid;
        using Source<TF>::fields;
};
#endif
