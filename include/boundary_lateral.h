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

#ifndef BOUNDARY_LATERAL_H
#define BOUNDARY_LATERAL_H

#include <vector>
#include <string>
#include <map>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;

template<typename TF>
class Boundary_lateral
{
    public:
        Boundary_lateral(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_lateral();

        void init();
        void create(Input&, const std::string&);
        void set_ghost_cells();
        void update_time_dependent(Timeloop<TF>&);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_inoutflow;
        bool sw_inoutflow_u;
        bool sw_inoutflow_v;
        std::vector<std::string> inoutflow_s;

        // Current (time interpolated) boundary conditions:
        std::map<std::string, std::vector<TF>> lbc_w;
        std::map<std::string, std::vector<TF>> lbc_s;
        std::map<std::string, std::vector<TF>> lbc_e;
        std::map<std::string, std::vector<TF>> lbc_n;

        // Input from NetCDF, for local MPI subdomain:
        std::vector<TF> time_in;
        std::map<std::string, std::vector<TF>> lbc_w_in;
        std::map<std::string, std::vector<TF>> lbc_s_in;
        std::map<std::string, std::vector<TF>> lbc_e_in;
        std::map<std::string, std::vector<TF>> lbc_n_in;
};
#endif
