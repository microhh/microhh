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

#ifndef SOIL_ENABLED
#define SOIL_ENABLED

#include "soil.h"
#include "soil_field.h"

//template<typename> class Soil_field;

template<typename TF>
class Soil_enabled : public Soil<TF>
{
    public:
        Soil_enabled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Soil_enabled();

        void init();
        void create_cold_start(Input&, Netcdf_handle&);
        void create_fields_grid_stats(Input&, Netcdf_handle&, Stats<TF>&, Cross<TF>&);
        void save_prognostic_fields(int);
        void load_prognostic_fields(int);
        void calc_tendencies() {};
        void exec_stats(Stats<TF>&);

    private:
        using Soil<TF>::sw_soil;
        using Soil<TF>::grid;
        using Soil<TF>::master;
        using Soil<TF>::fields;

        bool sw_interactive;
        bool sw_homogeneous;

        // Soil fields
        std::shared_ptr<Soil_field<TF>> t_soil;
        std::shared_ptr<Soil_field<TF>> theta_soil;

        // Soil properties
        std::vector<int> soil_index;    // Index in lookup tables

        std::vector<TF> diffusivity;    // Full level (m2 s-1)
        std::vector<TF> diffusivity_h;  // Half level (m2 s-1)
        std::vector<TF> conductivity;   // Full level (unit m s-1)
        std::vector<TF> conductivity_h; // Half level (unit m s-1)
        std::vector<TF> source;         // Source term (unit s-1)

        // Soil grid dimensions and data
        int ktot;
        std::vector<TF> z;
        std::vector<TF> zh;
        std::vector<TF> dz;
        std::vector<TF> dzh;
        std::vector<TF> dzi;
        std::vector<TF> dzhi;
        TF zsize;
};
#endif
