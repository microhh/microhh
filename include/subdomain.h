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

#ifndef SUBDOMAIN_H
#define SUBDOMAIN_H

#include <vector>
#include <string>
#include <map>

#include "master.h"
#include "grid.h"
#include "boundary_lateral_kernels.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;
template<typename> class Stats;
template<typename> class Field3d_io;


template<typename TF>
class Nn_data
{
    public:
        Nn_data() {}

        Nn_data(
            const std::vector<TF>& x_in,
            const std::vector<TF>& y_in,
            const std::vector<TF>& x,
            const std::vector<TF>& y,
            const Grid_data<TF>& gd,
            const MPI_data& md,
            const int ktot)
        {
            //this->itot_g = x_in.size();
            //this->jtot_g = y_in.size();
            //this->ktot_g = ktot;

            //// Bounds sub-domain on this MPI process.
            //const TF xsize_block = gd.imax * gd.dx;
            //const TF ysize_block = gd.jmax * gd.dy;

            //const TF xstart_block = md.mpicoordx * xsize_block;
            //const TF ystart_block = md.mpicoordy * ysize_block;

            //const TF xend_block = (md.mpicoordx + 1) * xsize_block;
            //const TF yend_block = (md.mpicoordy + 1) * ysize_block;

            //// Check if this MPI process contains data.
            //if (blk::intersects_mpi_subdomain(x_in, y_in, xstart_block, xend_block, ystart_block, yend_block))
            //{
            //    this->has_data = true;

            //    // Find range on current MPI process.
            //    this->i_range = blk::get_start_end_indexes(x_in, xstart_block, xend_block);
            //    this->j_range = blk::get_start_end_indexes(y_in, ystart_block, yend_block);

            //    // Slice out local coordinates, for finding NN indexes below.
            //    const std::vector<TF> x_in_s(x_in.begin() + this->i_range.first, x_in.begin() + this->i_range.second);
            //    const std::vector<TF> y_in_s(y_in.begin() + this->j_range.first, y_in.begin() + this->j_range.second);

            //    // Accounting.
            //    this->itot_s = i_range.second - i_range.first;
            //    this->jtot_s = j_range.second - j_range.first;
            //    this->ktot_s = this->ktot_g;

            //    this->istride = 1;
            //    this->jstride = this->itot_s;
            //    this->kstride = this->itot_s * this->jtot_s;

            //    // Find nearest-neighbour indexes.
            //    this->nn_i = blk::get_nn_indexes<TF>(x_in_s, x);
            //    this->nn_j = blk::get_nn_indexes<TF>(y_in_s, y);
            //}
            //else
            //{
            //    this->has_data = false;

            //    // Set sizes to zero for MPI-IO.
            //    this->i_range = std::make_pair<int>(0,0);
            //    this->j_range = std::make_pair<int>(0,0);

            //    this->itot_s = 0;
            //    this->jtot_s = 0;
            //    this->ktot_s = 0;
            //}

            //// Always resize, even if size is zero. This way we can get a valid pointer for MPI-IO.
            //this->fld.resize(itot_s * jtot_s * ktot_s);
        }

        // Not all tasks have data.
        bool has_data;

        // Global size (not accounting for MPI).
        int itot_g, jtot_g, ktot_g;

        // Start/end indices in global array. Needed for MPI-IO hyperslabs.
        std::pair<int, int> i_range;
        std::pair<int, int> j_range;

        // Local sub-set on this MPI process.
        int itot_s, jtot_s, ktot_s;
        int istride, jstride, kstride;

        // Data.
        std::vector<TF> fld;

        // Nearest neighbour indices.
        std::vector<int> nn_i;
        std::vector<int> nn_j;

        TF& operator()(const int i, const int j, const int k)
        {
            return fld[i*istride + j*jstride + k*kstride];
        }

        const TF& operator()(const int i, const int j, const int k) const
        {
            return fld[i*istride + j*jstride + k*kstride];
        }
};


template<typename TF>
class Subdomain
{
    public:
        Subdomain(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Subdomain();

        void create();
        unsigned long get_time_limit(unsigned long);
        void save_bcs(Timeloop<TF>&);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_io<TF> field3d_io;

        bool sw_subdomain;
        bool sw_save_wtop;
        bool sw_save_buffer;

        TF xstart_sub, xend_sub;
        TF ystart_sub, yend_sub;
        TF xsize_sub, ysize_sub;
        TF dx_sub, dy_sub;

        int itot_sub, jtot_sub;

        int grid_ratio_sub;
        int n_ghost_sub;
        int n_sponge_sub;

        // NOTE: savetime_sub is always used for both LBCs and w_top
        //       These NEED to be on the same frequency.
        //       The 3D sponge layer is optionally saved at a lower frequency.
        unsigned int savetime_bcs;
        unsigned int savetime_buffer;

        std::map<std::string, Nn_data<TF>> lbc_sub_w;
        std::map<std::string, Nn_data<TF>> lbc_sub_e;
        std::map<std::string, Nn_data<TF>> lbc_sub_s;
        std::map<std::string, Nn_data<TF>> lbc_sub_n;

        Nn_data<TF> bc_wtop_sub;

        TF zstart_buffer;
        int buffer_kstart;
        int buffer_kstarth;

        Nn_data<TF> bc_buffer;
        Nn_data<TF> bc_bufferh;
};
#endif