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

#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <vector>
#include <string>
#include <map>

#include "timeloop.h"
#include "boundary.h"
#include "stats.h"

enum Jval
{
    o31d   = 0,
    h2o2   = 1,
    no2    = 2,
    no3    = 3,
    n2o5   = 4,
    ch2or  = 5,
    ch2om  = 6,
    ch3o2h = 7,
    n_jval
};

class Master;
class Input;
template<typename> class Grid;
template<typename> class Soil_grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Deposition;

/**
 * Class that creates a chemistry term for scalars.
 */

enum class Chemistry_type {disabled, enabled, simple};

template<typename TF>
class Chemistry
{
    public:
        Chemistry(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the chemistry class.
        ~Chemistry();                                       ///< Destructor  of the chemistry class.

        void init(Input&);                                                                      ///< Initialize the arrays that contain the profiles.
        void create(const Timeloop<TF>&, std::string, Netcdf_handle&, Stats<TF>&, Cross<TF>&);  ///< Read the profiles of the forces from the input.
        void update_time_dependent(Timeloop<TF>&,Boundary<TF>&);                                ///< Update the time dependent parameters.
        void exec(Thermo<TF>&, double, double);                                                 ///< Add the tendencies belonging to the chemistry processes.
        void exec_stats(const int, const double, Stats<TF>&);                                   ///< Execute statistics
        void exec_cross(Cross<TF>&, unsigned long);                                             ///< Execute cross-sections.

    protected:
        // Cross sections
        std::vector<std::string> cross_list;         // List of active cross variables

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_chemistry;

        Field3d_operators<TF> field3d_operators;

        std::shared_ptr<Deposition<TF>> deposition;

        std::vector<std::string> jname={"jo31d","jh2o2","jno2","jno3","jn2o5","jch2or","jch2om","jch3o2h"};

        std::vector<TF> jval;   // time-interpolated value to pass to the chemistry routine

        // CPU arrays.
        std::vector<double> time;   // NOTE: keep this double.
        std::vector<TF> jo31d;
        std::vector<TF> jh2o2;
        std::vector<TF> jno2;
        std::vector<TF> jno3;
        std::vector<TF> jn2o5;
        std::vector<TF> jch2or;
        std::vector<TF> jch2om;
        std::vector<TF> jch3o2h;

        // vectors to contain calculated deposition velocities (m/s)
        std::vector<TF> vdo3;
        std::vector<TF> vdno;
        std::vector<TF> vdno2;
        std::vector<TF> vdhno3;
        std::vector<TF> vdh2o2;
        std::vector<TF> vdrooh;
        std::vector<TF> vdhcho;


};
#endif
