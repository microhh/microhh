/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef FIELD3D_H
#define FIELD3D_H

#include <string>
#include <vector>
#include <array>

#include "cuda_buffer.h"

class Master;
template<typename> class Grid;

template<typename TF>
class Field3d
{
    public:
        // Functions
        Field3d(
                Master&, Grid<TF>&,
                const std::string&, const std::string&, const std::string&, const std::string&,
                const std::array<int,3>&);
        ~Field3d();

        int init();

        // Variables at CPU.
        std::vector<TF> fld;
        std::vector<TF> fld_bot;
        std::vector<TF> fld_top;
        std::vector<TF> fld_mean;
        std::vector<TF> grad_bot;
        std::vector<TF> grad_top;
        std::vector<TF> flux_bot;
        std::vector<TF> flux_top;

        std::string name;
        std::string unit;
        std::string longname;
        std::string group;

        std::array<int,3> loc;

        TF visc;

        // Device functions and variables
        void init_device();  // Allocate Field3D fields at device
        void clear_device(); // Deallocate Field3D fields at device

        cuda_vector<TF> fld_g;
        cuda_vector<TF> fld_bot_g;
        cuda_vector<TF> fld_top_g;
        cuda_vector<TF> fld_mean_g;
        cuda_vector<TF> grad_bot_g;
        cuda_vector<TF> grad_top_g;
        cuda_vector<TF> flux_bot_g;
        cuda_vector<TF> flux_top_g;

    private:
        Master& master;
        Grid<TF>& grid;
};
#endif
