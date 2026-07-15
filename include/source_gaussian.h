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

#ifndef SOURCE_GAUSSIAN_H
#define SOURCE_GAUSSIAN_H

#include <set>

#include "source.h"

class Input;

template<typename TF>
class Source_gaussian : public Source<TF>
{
    public:
        Source_gaussian(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Source_gaussian();

        void init();
        void create(Input&, Timeloop<TF>&, Netcdf_handle&);
        void exec(Thermo<TF>&, Timeloop<TF>&);
        void update_time_dependent(Timeloop<TF>&);

        #ifdef USECUDA
        void prepare_device();
        #endif

    private:
        using Source<TF>::master;
        using Source<TF>::grid;
        using Source<TF>::fields;

        struct Shape
        {
            std::vector<int> range_x;
            std::vector<int> range_y;
            std::vector<int> range_z;
        };

        std::vector<Shape> shape;

        bool swsource;
        std::vector<bool> sw_vmr;

        TF x0;
        TF y0;
        TF z0;

        std::vector<std::string> sourcelist;
        std::vector<TF> source_x0;
        std::vector<TF> source_y0;
        std::vector<TF> source_z0;
        std::vector<TF> sigma_x;
        std::vector<TF> sigma_y;
        std::vector<TF> sigma_z;
        std::vector<TF> strength;
        std::vector<TF> line_x;
        std::vector<TF> line_y;
        std::vector<TF> line_z;
        std::vector<TF> norm;

        // Timedep source location and strength
        bool swtimedep_location;
        bool swtimedep_strength;
        std::map<std::string, Timedep<TF>*> tdep_source_x0;
        std::map<std::string, Timedep<TF>*> tdep_source_y0;
        std::map<std::string, Timedep<TF>*> tdep_source_z0;
        std::map<std::string, Timedep<TF>*> tdep_source_strength;

        // Option for vertical profiles
        bool sw_emission_profile;
        std::vector<int> profile_index;
        std::set<int> unique_profile_indexes;
        std::map<int, std::vector<TF>> profile_z;

        TF calc_norm(
                const TF* const, const TF, const TF, const TF,
                const TF* const, const TF, const TF, const TF,
                const TF* const, const TF, const TF, const TF,
                const TF* const,
                std::vector<int>, std::vector<int>, std::vector<int>,
                const TF* const, const bool, const bool);
};
#endif
