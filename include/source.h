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

#ifndef SOURCE_H
#define SOURCE_H

#include <memory>

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;


template<typename TF>
class Source
{
    public:
        Source(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Source();

        static std::shared_ptr<Source> factory(Master&, Grid<TF>&, Fields<TF>&, Input&);

        virtual void init() = 0;
        virtual void create(Input&, Timeloop<TF>&, Netcdf_handle&) = 0;
        virtual void exec() = 0;
        virtual void update_time_dependent(Timeloop<TF>&) = 0;

        #ifdef USECUDA
        virtual void prepare_device() = 0;
        #endif

    protected:
        Master& master; ///< Pointer to master class.
        Grid<TF>& grid; ///< Pointer to grid class.
        Fields<TF>& fields; ///< Pointer to fields class.
};
#endif
