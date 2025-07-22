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

#ifndef SUBDOMAINS_H
#define SUBDOMAINS_H

#include <vector>
#include <string>
#include <map>

#include "master.h"
#include "grid.h"
#include "nn_interpolator.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;
template<typename> class Stats;
template<typename> class Field3d_io;


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

        TF xstart, xend;
        TF ystart, yend;
        TF xsize, ysize;
        TF dx, dy;

        int itot, jtot;

        int grid_ratio_ij;
        int grid_ratio_k;

        int n_ghost;
        int n_sponge;

        // NOTE: savetime is always used for both LBCs and w_top
        //       These NEED to be on the same frequency.
        //       The 3D sponge layer is optionally saved at a different frequency.
        unsigned int savetime_bcs;
        unsigned int savetime_buffer;

        std::map<std::string, std::unique_ptr<NN_interpolator<TF>>> lbc_w;
        std::map<std::string, std::unique_ptr<NN_interpolator<TF>>> lbc_e;
        std::map<std::string, std::unique_ptr<NN_interpolator<TF>>> lbc_s;
        std::map<std::string, std::unique_ptr<NN_interpolator<TF>>> lbc_n;

        std::unique_ptr<NN_interpolator<TF>> bc_wtop;

        TF zstart_buffer;
        int buffer_kstart;
        int buffer_kstarth;

        std::unique_ptr<NN_interpolator<TF>> bc_buffer;
        std::unique_ptr<NN_interpolator<TF>> bc_bufferh;
};
#endif
