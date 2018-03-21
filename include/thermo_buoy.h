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

#ifndef THERMO_BUOY
#define THERMO_BUOY

#include "thermo.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
class Data_block;



/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */
template<typename TF>
class Thermo_buoy : public Thermo<TF>
{
    public:
        Thermo_buoy(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the dry thermodynamics class.
        virtual ~Thermo_buoy();        ///< Destructor of the dry thermodynamics class.

        void exec(); ///< Add the tendencies belonging to the buoyancy.
        unsigned long get_time_limit(unsigned long, double); ///< Compute the time limit (n/a for thermo_buoy)

        bool check_field_exists(std::string name);
        void get_buoyancy_surf(Field3d<TF>&);             ///< Compute the near-surface and bottom buoyancy for usage in another routine.
        void get_buoyancy_fluxbot(Field3d<TF>&);           ///< Compute the bottom buoyancy flux for usage in another routine.
        void get_prog_vars(std::vector<std::string>&); ///< Retrieve a list of prognostic variables.
        void get_thermo_field(Field3d<TF>&, std::string, bool); ///< Compute the buoyancy for usage in another routine.
        TF get_buoyancy_diffusivity();

        // Empty functions that are allowed to pass.
        void init() {}
        void create(Input&, Data_block&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&) {}
        void exec_stats(Stats<TF>&, std::string, Field3d<TF>&, Field3d<TF>&, const Diff<TF>&) {};
        void exec_cross(Cross<TF>&, unsigned long) {};
        void exec_dump(Dump<TF>&, unsigned long) {};
        void exec_column(Column<TF>&) {};
        void get_mask(Field3d<TF>&, Field3d<TF>&, Stats<TF>&, std::string) {};
        void update_time_dependent() {}

#ifdef USECUDA
        void prepare_device() {};
        void clear_device() {};
        void forward_device() {};
        void backward_device() {};
#endif

private:
        using Thermo<TF>::swthermo;
        using Thermo<TF>::master;
        using Thermo<TF>::grid;
        using Thermo<TF>::fields;
        struct background_state
        {
            TF alpha;  ///< Slope angle in radians.
            TF n2;     ///< Background stratification.
            bool has_slope; ///< Boolean switch for slope flows
            bool has_N2;    ///< Boolean switch for imposed stratification
        };
        background_state bs;
        background_state bs_stats;
};
#endif
