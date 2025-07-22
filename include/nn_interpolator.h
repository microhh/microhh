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

#ifndef NN_INTERPOLATE_H
#define NN_INTERPOLATE_H

#include <vector>

#include "master.h"
#include "grid.h"
#include "nn_interpolator_kernels.h"

namespace nn = NN_interpolator_kernels;


template<typename TF>
class NN_interpolator
{
    public:
        NN_interpolator() {}

        NN_interpolator(
            const std::vector<TF>& x_in,
            const std::vector<TF>& y_in,
            const std::vector<TF>& z_in,
            const std::vector<TF>& x,
            const std::vector<TF>& y,
            const std::vector<TF>& z,
            const Grid_data<TF>& gd,
            const MPI_data& md);

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
        std::vector<int> nn_k;

        TF& operator()(const int i, const int j, const int k)
        {
            return fld[i*istride + j*jstride + k*kstride];
        }

        const TF& operator()(const int i, const int j, const int k) const
        {
            return fld[i*istride + j*jstride + k*kstride];
        }
};
#endif
