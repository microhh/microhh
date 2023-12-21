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
template<typename> class Field3d_io;

enum class Lbc_location {West, East, South, North};

template<typename TF>
class Boundary_lateral
{
    public:
        Boundary_lateral(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_lateral();

        void init();
        void create(Input&, Timeloop<TF>&, const std::string&);
        void set_ghost_cells(Timeloop<TF>&);
        void update_time_dependent(Timeloop<TF>&, const bool pres_fix=false);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_io<TF> field3d_io;

        void read_xy_slice(
                std::vector<TF>&, const std::string&, const int);

        bool sw_inoutflow;
        bool sw_inoutflow_uv;
        bool sw_inoutflow_w;
        bool sw_neumann_w;
        bool sw_timedep;
        bool sw_wtop_2d;

        std::vector<std::string> inoutflow_s;

        // Sponge/diffusion layer:
        bool sw_sponge;
        int n_sponge;
        TF tau_nudge;
        TF w_diff;

        // Boundary perturbations:
        bool sw_perturb;
        int perturb_width;
        int perturb_block;
        int perturb_seed;
        int perturb_kend;
        std::vector<std::string> perturb_list;
        std::map<std::string, TF> perturb_ampl;

        // Turbulence recycling.
        bool sw_recycle;
        std::vector<std::string> recycle_list;
        TF tau_recycle;
        int recycle_offset;

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

        std::vector<TF> div_u;
        std::vector<TF> div_v;
        std::vector<TF> w_top_in;
        std::vector<TF> w_top;

        unsigned int lbc_load_freq;
        unsigned long itime_w_top_prev;
        unsigned long itime_w_top_next;

        std::vector<TF> w_top_prev;
        std::vector<TF> w_top_next;
};
#endif
