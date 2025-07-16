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

#ifndef BOUNDARY_LATERAL_KERNELS_H
#define BOUNDARY_LATERAL_KERNELS_H

#include <vector>


namespace boundary_lateral_kernels
{
    template<typename TF>
    std::vector<TF> arange(
        const TF start,
        const TF stop,
        const TF step)
    {
        // Determine size.
        int size;
        if (step > 0 && start < stop)
            size = static_cast<int>((stop - start + step - 1) / step);
        else if (step < 0 && start > stop)
            size = static_cast<int>((start - stop - step - 1) / (-step));
        else
            throw std::runtime_error("Invalid configuration for arange()!");

        // Allocate vector.
        std::vector<TF> result;
        result.reserve(size);

        // Fill vector.
        if (step > 0)
            for (TF value = start; value < stop; value += step)
                result.push_back(value);
        else
            for (TF value = start; value > stop; value += step)
                result.push_back(value);

        return result;
    }


    template<typename TF>
    bool intersects_mpi_subdomain(
        const std::vector<TF>& x,
        const std::vector<TF>& y,
        const TF x_min, const TF x_max,
        const TF y_min, const TF y_max)
    {
        // First simple check on bounds.
        if (x.back() < x_min || x.front() > x_max)
            return false;

        if (y.back() < y_min || y.front() > y_max)
            return false;

        // Check coordinate pairs.
        for (int i=0; i<x.size(); ++i)
            for (int j=0; j<y.size(); ++j)
                if (x[i] >= x_min && x[i] < x_max && y[j] >= y_min && y[j] < y_max)
                    return true;

        return false;
    }


    template<typename TF>
    std::pair<int, int> get_start_end_indexes(
        const std::vector<TF>& x,
        const TF x0,
        const TF x1)
    {
        auto start_it = std::lower_bound(x.begin(), x.end(), x0);
        auto end_it = std::lower_bound(x.begin(), x.end(), x1);

        int start_idx = std::distance(x.begin(), start_it);
        int end_idx = std::distance(x.begin(), end_it);

        return std::make_pair(start_idx, end_idx);
    }


    template<typename TF>
    std::vector<int> get_nn_indexes(
        const std::vector<TF>& x_target,
        const std::vector<TF>& x_all)
    {
        std::vector<int> nn_index(x_target.size());
        int nearest_index;

        for (int i=0; i<nn_index.size(); ++i)
        {
            auto it = std::lower_bound(x_all.begin(), x_all.end(), x_target[i]);

            if (it == x_all.begin())
                nearest_index = 0;
            else if (it == x_all.end())
                nearest_index = x_all.size() - 1;
            else
            {
                // std::lower_bound returns the first element that is larger, not the nearest...
                const int upper_index = std::distance(x_all.begin(), it);
                const int lower_index = upper_index - 1;

                TF dist_to_lower = x_target[i] - x_all[lower_index];
                TF dist_to_upper = x_all[upper_index] - x_target[i];

                nearest_index = (dist_to_lower <= dist_to_upper) ? lower_index : upper_index;
            }

            nn_index[i] = nearest_index;
        }

        return nn_index;
    }


    template<typename TF>
    void nn_interpolate(
        TF* const restrict fld_out,
        const TF* const restrict fld_in,
        const int* const restrict nn_i,
        const int* const restrict nn_j,
        const int isize,
        const int jsize,
        const int ksize,
        const int kgc,
        const int istride_out,
        const int jstride_out,
        const int kstride_out,
        const int istride_in,
        const int jstride_in,
        const int kstride_in)
    {
        for (int k=0; k<ksize; ++k)
            for (int j=0; j<jsize; ++j)
                for (int i=0; i<isize; ++i)
                {
                    const int ijk_in = nn_i[i]*istride_in + nn_j[j]*jstride_in + (k+kgc)*kstride_in;
                    const int ijk_out = i*istride_out + j*jstride_out + k*kstride_out;

                    fld_out[ijk_out] = fld_in[ijk_in];
                }
    }


}
#endif
