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


// Class that holds a single LBC.
template<typename TF>
class Lbc
{
    public:
        Lbc() {}

        Lbc(
                const int dim_x,
                const int dim_y,
                const int dim_z)
        {
            this->dim_x = dim_x;
            this->dim_y = dim_y;
            this->dim_z = dim_z;

            size = dim_x * dim_y * dim_z;

            x_stride = 1;
            y_stride = dim_x;
            z_stride = dim_x * dim_y;

            data = std::vector<TF>(size, TF(0));
        }

        int dim_x, dim_y, dim_z;
        int size;
        int x_stride, y_stride, z_stride;

        std::vector<TF> data;

        TF& operator()(const int i, const int j, const int k)
        {
            return data[i * x_stride + j * y_stride + k * z_stride];
        }

        const TF& operator()(const int i, const int j, const int k) const
        {
            return data[i * x_stride + j * y_stride + k * z_stride];
        }
};


/* 
 * Class that holds all the lateral boundary conditions.
 * Wrapped in a class, since we can have multiple instances of this class (own domain, child, etc.).
*/
template<typename TF>
class Lbcs
{
    public:
        Lbcs() {}

        Lbcs(
                const std::vector<std::string>& scalar_list,
                const int itot, const int jtot, const int ktot,
                const int n_ghost, const int n_sponge)
        {
            this->itot = itot;
            this->jtot = jtot;
            this->ktot = ktot;

            this->n_ghost = n_ghost;
            this->n_sponge = n_sponge;
            n_lbc = n_ghost + n_sponge;

            // Dimensions of all possible LBCs.
            dim_x = itot + 2*n_ghost;
            dim_xgw = itot + 2*n_ghost;
            dim_xge = n_lbc;
            dim_xhgw = n_lbc + 1;
            dim_xhge = n_lbc;

            dim_y = jtot + 2*n_ghost;
            dim_yh = jtot + 2*n_ghost;
            dim_ygs = n_lbc;
            dim_ygn = n_lbc;
            dim_yhgs = n_lbc + 1;
            dim_yhgn = n_lbc;

            dim_z = ktot;
            dim_zh = ktot;

            // Momentum fields.
            lbc_w["u"] = Lbc<TF>(dim_xhgw, dim_y,   dim_z);
            lbc_e["u"] = Lbc<TF>(dim_xhge, dim_y,   dim_z);
            lbc_s["u"] = Lbc<TF>(dim_xh,   dim_ygs, dim_z);
            lbc_n["u"] = Lbc<TF>(dim_xh,   dim_ygn, dim_z);

            lbc_w["v"] = Lbc<TF>(dim_xgw, dim_yh,   dim_z);
            lbc_e["v"] = Lbc<TF>(dim_xge, dim_yh,   dim_z);
            lbc_s["v"] = Lbc<TF>(dim_x,   dim_yhgs, dim_z);
            lbc_n["v"] = Lbc<TF>(dim_x,   dim_yhgn, dim_z);

            lbc_w["w"] = Lbc<TF>(dim_xgw, dim_y,   dim_zh);
            lbc_e["w"] = Lbc<TF>(dim_xge, dim_y,   dim_zh);
            lbc_s["w"] = Lbc<TF>(dim_x,   dim_ygs, dim_zh);
            lbc_n["w"] = Lbc<TF>(dim_x,   dim_ygn, dim_zh);

            // Scalar fields.
            for (const auto& scalar : scalar_list)
            {
                lbc_w[scalar] = Lbc<TF>(dim_xgw, dim_y,   dim_z);
                lbc_e[scalar] = Lbc<TF>(dim_xge, dim_y,   dim_z);
                lbc_s[scalar] = Lbc<TF>(dim_x,   dim_ygs, dim_z);
                lbc_n[scalar] = Lbc<TF>(dim_x,   dim_ygn, dim_z);
            }
        }

        // `Lbc` instances at all domain edges.
        std::map<std::string, Lbc<TF>> lbc_n;
        std::map<std::string, Lbc<TF>> lbc_e;
        std::map<std::string, Lbc<TF>> lbc_s;
        std::map<std::string, Lbc<TF>> lbc_w;

        int itot, jtot, ktot;
        int n_ghost, n_sponge, n_lbc;

        // Dimension sizes.
        int dim_x, dim_xh, dim_xgw, dim_xge, dim_xhgw, dim_xhge;
        int dim_y, dim_yh, dim_ygs, dim_ygn, dim_yhgs, dim_yhgn;
        int dim_z, dim_zh;
};


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
        std::map<Lbc_location, bool> sw_recycle;
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

        Lbcs<TF> lbcs_sub;
};
#endif