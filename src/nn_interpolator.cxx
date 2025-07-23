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

#include <vector>
#include <algorithm>
#include <cmath>

#include "nn_interpolator.h"
#include "nn_interpolator_kernels.h"

#include "master.h"
#include "grid.h"

namespace nnk = NN_interpolator_kernels;

namespace
{
    #ifdef USEMPI
    template<typename TF> MPI_Datatype mpi_fp_type();
    template<> MPI_Datatype mpi_fp_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_fp_type<float>() { return MPI_FLOAT; }
    #endif
}


template<typename TF>
NN_interpolator<TF>::NN_interpolator(
    const std::vector<TF>& x_in,
    const std::vector<TF>& y_in,
    const std::vector<TF>& z_in,
    const std::vector<TF>& x,
    const std::vector<TF>& y,
    const std::vector<TF>& z,
    Master& masterin,
    Grid<TF>& gridin) : master(masterin), grid(gridin)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    this->itot_g = x_in.size();
    this->jtot_g = y_in.size();
    this->ktot_g = z_in.size();

    // Bounds sub-domain on this MPI process.
    const TF xsize_block = gd.imax * gd.dx;
    const TF ysize_block = gd.jmax * gd.dy;

    const TF xstart_block = md.mpicoordx * xsize_block;
    const TF ystart_block = md.mpicoordy * ysize_block;

    const TF xend_block = (md.mpicoordx + 1) * xsize_block;
    const TF yend_block = (md.mpicoordy + 1) * ysize_block;

    // Check if this MPI process contains data.
    if (nnk::intersects_mpi_subdomain(x_in, y_in, xstart_block, xend_block, ystart_block, yend_block))
    {
        this->has_data = true;

        // Find range on current MPI process.
        this->i_range = nnk::get_start_end_indexes(x_in, xstart_block, xend_block);
        this->j_range = nnk::get_start_end_indexes(y_in, ystart_block, yend_block);

        // Slice out local coordinates, for finding NN indexes below.
        const std::vector<TF> x_in_s(x_in.begin() + this->i_range.first, x_in.begin() + this->i_range.second);
        const std::vector<TF> y_in_s(y_in.begin() + this->j_range.first, y_in.begin() + this->j_range.second);

        // Accounting.
        this->itot_s = i_range.second - i_range.first;
        this->jtot_s = j_range.second - j_range.first;
        this->ktot_s = this->ktot_g;

        this->istride = 1;
        this->jstride = this->itot_s;
        this->kstride = this->itot_s * this->jtot_s;

        // Find nearest-neighbour indexes.
        this->nn_i = nnk::get_nn_indexes<TF>(x_in_s, x);
        this->nn_j = nnk::get_nn_indexes<TF>(y_in_s, y);
        this->nn_k = nnk::get_nn_indexes<TF>(z_in,   z);

        //std::cout << "i: range=" << i_range.first << ":" << i_range.second << std::endl;

        //std::cout << "x: ";
        //for (auto& x : x_in)
        //    std::cout << x << ", ";
        //std::cout << std::endl;

        //std::cout << "nn_i: ";
        //for (auto& i : nn_i)
        //    std::cout << i << ", ";
        //std::cout << std::endl;

        //std::cout << "j: range=" << j_range.first << ":" << j_range.second << std::endl;

        //std::cout << "y: ";
        //for (auto& y : y_in)
        //    std::cout << y << ", ";
        //std::cout << std::endl;

        //std::cout << "nn_j: ";
        //for (auto& j : nn_j)
        //    std::cout << j << ", ";
        //std::cout << std::endl;

        //std::cout << "k:" << std::endl;

        //std::cout << "z: ";
        //for (auto& z : z_in)
        //    std::cout << z << ", ";
        //std::cout << std::endl;

        //std::cout << "nn_k: ";
        //for (auto& k : nn_k)
        //    std::cout << k << ", ";
        //std::cout << std::endl;

        //std::cout << "z:" << std::endl;
        //for (int k=0; k<z.size(); ++k)
        //    std::cout << k << ": " << z[k] << std::endl;

        //throw 1;
    }
    else
    {
        this->has_data = false;

        // Set sizes to zero for MPI-IO.
        this->i_range = std::make_pair<int>(0,0);
        this->j_range = std::make_pair<int>(0,0);

        this->itot_s = 0;
        this->jtot_s = 0;
        this->ktot_s = 0;
    }

    // Always resize, even if size is zero. This way we can get a valid pointer for MPI-IO.
    this->fld.resize(itot_s * jtot_s * ktot_s);
}


template<typename TF>
void NN_interpolator<TF>::interpolate(const std::vector<TF>& fld_in)
{
    auto& gd = grid.get_grid_data();

    nnk::nn_interpolate(
        this->fld.data(),
        fld_in.data(),
        this->nn_i.data(),
        this->nn_j.data(),
        this->nn_k.data(),
        this->itot_s,
        this->jtot_s,
        this->ktot_s,
        this->istride,
        this->jstride,
        this->kstride,
        gd.istride,
        gd.jstride,
        gd.kstride);
}


#ifdef FLOAT_SINGLE
template class NN_interpolator<float>;
#else
template class NN_interpolator<double>;
#endif
