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

#ifndef SOURCE_3D_H
#define SOURCE_3D_H

#include "source.h"

class Input;

template<typename TF>
class Source_3d : public Source<TF>
{
    public:
        Source_3d(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Source_3d();

        void init();
        void create(Input&, Netcdf_handle&);
        void exec();
        void update_time_dependent(Timeloop<TF>&);

    private:
        using Source<TF>::master;
        using Source<TF>::grid;
        using Source<TF>::fields;

        std::vector<std::string> sourcelist;
        std::map<std::string, std::vector<TF>> emission;

        bool sw_timedep;    // Switch for time dependent 3D input.
        int ktot;           // Number of vertical levels with emissions.
};
#endif