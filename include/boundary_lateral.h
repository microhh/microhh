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
template<typename> class Stats;
template<typename> class Field3d_io;

enum class Lbc_location {West, East, South, North};

template<typename TF>
using Lbc_map = std::map<std::string, std::vector<TF>>;

template<typename TF>
class Boundary_lateral
{
    public:
        Boundary_lateral(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_lateral();

        void init();
        void create(Input&, Timeloop<TF>&, Stats<TF>&, const std::string&);
        void set_ghost_cells(Timeloop<TF>&);
        void exec_lateral_sponge(Stats<TF>&);
        void update_time_dependent(Timeloop<TF>&, const bool pres_fix=false);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_io<TF> field3d_io;

        void read_lbc(TF&, TF&, Lbc_map<TF>&, Lbc_map<TF>&, Lbc_map<TF>&, Lbc_map<TF>&, const int);
        void read_xy_slice(
                std::vector<TF>&, const std::string&, const int);

        bool sw_openbc;
        bool sw_openbc_uv;
        bool sw_openbc_w;
        bool sw_neumann_w;
        bool sw_timedep;
        bool sw_wtop_2d;

        std::vector<std::string> slist;

        // Sponge/diffusion layer:
        bool sw_sponge;
        int n_sponge;
        TF tau_sponge;
        TF w_diff;

        // Turbulence recycling.
        bool sw_recycle;
        std::vector<std::string> recycle_list;
        TF tau_recycle;
        int recycle_offset;

        // Current (constant or time interpolated) BCs:
        Lbc_map<TF> lbc_w;
        Lbc_map<TF> lbc_s;
        Lbc_map<TF> lbc_e;
        Lbc_map<TF> lbc_n;

        // Previous and next BCs, for time dependent LBCs.
        Lbc_map<TF> lbc_w_prev;
        Lbc_map<TF> lbc_e_prev;
        Lbc_map<TF> lbc_s_prev;
        Lbc_map<TF> lbc_n_prev;

        Lbc_map<TF> lbc_w_next;
        Lbc_map<TF> lbc_e_next;
        Lbc_map<TF> lbc_s_next;
        Lbc_map<TF> lbc_n_next;

        unsigned int loadfreq;
        unsigned long prev_itime;
        unsigned long next_itime;

        // Spatially constant (calculated) `w_top`.
        TF w_top;
        TF w_top_prev;
        TF w_top_next;

        // Spatially varying (input) `w_top`.
        std::vector<TF> w_top_2d;
        std::vector<TF> w_top_2d_prev;
        std::vector<TF> w_top_2d_next;

        const std::string tend_name = "lbc_sponge";
        const std::string tend_longname = "Lateral sponge layer";

        // Lateral boundary output for child domain.
        bool sw_subdomain;
        TF xstart_sub;
        TF xend_sub;
        TF ystart_sub;
        TF yend_sub;
        int refinement_sub;
        int n_ghost_sub;
        int n_sponge_sub;
};
#endif