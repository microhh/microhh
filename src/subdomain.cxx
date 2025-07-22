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

#include "subdomain.h"

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "constants.h"
#include "nn_interpolator_kernels.h"

namespace nnk = NN_interpolator_kernels;


namespace
{
    #ifdef USEMPI
    template<typename TF> MPI_Datatype mpi_fp_type();
    template<> MPI_Datatype mpi_fp_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_fp_type<float>() { return MPI_FLOAT; }
    #endif


    template<typename TF>
    bool is_equal(const TF a, const TF b)
    {
        const TF epsilon = std::max(TF(10) * std::numeric_limits<TF>::epsilon(), std::max(std::abs(a), std::abs(b)) * TF(10) * std::numeric_limits<TF>::epsilon());
        return std::abs(a - b) <= epsilon;
    }
}


template<typename TF>
Subdomain<TF>::Subdomain(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin), field3d_io(masterin, gridin)
{
    sw_subdomain = inputin.get_item<bool>("subdomain", "sw_subdomain", "", false);

    if (sw_subdomain)
    {
        sw_save_wtop = inputin.get_item<bool>("subdomain", "sw_save_wtop", "", false);
        sw_save_buffer = inputin.get_item<bool>("subdomain", "sw_save_buffer", "", false);

        xstart = inputin.get_item<TF>("subdomain", "xstart", "");
        xend   = inputin.get_item<TF>("subdomain", "xend", "");
        ystart = inputin.get_item<TF>("subdomain", "ystart", "");
        yend   = inputin.get_item<TF>("subdomain", "yend", "");

        grid_ratio_ij = inputin.get_item<int>("subdomain", "grid_ratio_ij", "");
        grid_ratio_k = inputin.get_item<int>("subdomain", "grid_ratio_k", "", 1);

        n_ghost = inputin.get_item<int>("subdomain", "n_ghost", "");
        n_sponge = inputin.get_item<int>("subdomain", "n_sponge", "");

        savetime_bcs = inputin.get_item<int>("subdomain", "savetime_bcs", "");

        if (sw_save_buffer)
        {
            savetime_buffer = inputin.get_item<int>("subdomain", "savetime_buffer", "");
            zstart_buffer = inputin.get_item<TF>("subdomain", "zstart_buffer", "");
        }

        // Checks.
        if (grid_ratio_ij < 1 || grid_ratio_k < 1)
            throw std::runtime_error("Grid refinement ratios should be equal to or larger than 1.");
    }
}


template <typename TF>
Subdomain<TF>::~Subdomain()
{
}


template <typename TF>
void Subdomain<TF>::create()
{
    if (!sw_subdomain)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    bool error = false;

    auto check_subdomain = [&](const TF value, const TF spacing, const std::string& name)
    {
        if (!is_equal(std::fmod(value, spacing), TF(0)))
        {
            error = true;
            master.print_message(
                "ERROR: " + name + " must be an integer multiple of " + std::to_string(spacing) + ".");
        }
    };

    // Sub-domain has perfectly align with parent grid at x/y half levels.
    check_subdomain(xstart, gd.dx, "xstart");
    check_subdomain(xend, gd.dx, "xend");
    check_subdomain(ystart, gd.dy, "ystart");
    check_subdomain(yend, gd.dy, "yend");

    if (error)
        throw std::runtime_error("Sub-domain boundaries not aligned with parent grid.");

    // Initialise LBCs instance.
    xsize = xend - xstart;
    ysize = yend - ystart;

    itot = static_cast<int>(xsize / gd.dx * grid_ratio_ij + 0.5);
    jtot = static_cast<int>(ysize / gd.dy * grid_ratio_ij + 0.5);

    dx = xsize / itot;
    dy = ysize / jtot;

    // Coordinates of LBCs at different grid locations.
    // Full levels x/y.
    std::vector<TF> x_ns = nnk::arange<TF>(xstart - (n_ghost  - 0.5) * dx, xend   + n_ghost  * dx, dx);
    std::vector<TF> x_w  = nnk::arange<TF>(xstart - (n_ghost  - 0.5) * dx, xstart + n_sponge * dx, dx);
    std::vector<TF> x_e  = nnk::arange<TF>(xend   - (n_sponge - 0.5) * dx, xend   + n_ghost  * dx, dx);

    std::vector<TF> y_ew = nnk::arange<TF>(ystart - (n_ghost  - 0.5) * dy, yend   + n_ghost  * dy, dy);
    std::vector<TF> y_n  = nnk::arange<TF>(yend   - (n_sponge - 0.5) * dy, yend   + n_ghost  * dy, dy);
    std::vector<TF> y_s  = nnk::arange<TF>(ystart - (n_ghost  - 0.5) * dy, ystart + n_sponge * dy, dy);

    // Half level `x` for `u`.
    std::vector<TF> xh_ns = nnk::arange<TF>(xstart - n_ghost  * dx, xend   + n_ghost  * dx,        dx);
    std::vector<TF> xh_w  = nnk::arange<TF>(xstart - n_ghost  * dx, xstart + n_sponge * dx + dx/2, dx);
    std::vector<TF> xh_e  = nnk::arange<TF>(xend   - n_sponge * dx, xend   + n_ghost  * dx,        dx);

    // Half level `y` for `v`.
    std::vector<TF> yh_ew = nnk::arange<TF>(ystart - n_ghost  * dy, yend   + n_ghost  * dy,        dy);
    std::vector<TF> yh_n  = nnk::arange<TF>(yend   - n_sponge * dy, yend   + n_ghost  * dy,        dy);
    std::vector<TF> yh_s  = nnk::arange<TF>(ystart - n_ghost  * dy, ystart + n_sponge * dy + dy/2, dy);

    // Vertical coordinates are more tricky; by refining the grid, we cut each full level confined by
    // `zh[k] -> zh[k+1]` into `grid_ratio_k` equal layers, and next calculate the full level heights
    // that satisfy our grid definition `z[k] = 0.5 * (zh[k-1] + zh[k])`. For strongly stretched grids,
    // this results in strange locations for `z`...
    std::vector<TF> z;
    std::vector<TF> zh;

    if (grid_ratio_k > 1)
    {
        const int ktot_n = grid_ratio_k * gd.ktot;

        // Define new grid location (without ghost cells).
        z.resize(ktot_n);
        std::vector<TF> zh_tmp(ktot_n+1);

        // Calculate new half level heights.
        for (int k=0; k<gd.ktot; ++k)
        {
            const int kk = k + gd.kgc;
            const TF dz = gd.dz[kk] / grid_ratio_k;

            for (int s=0; s<grid_ratio_k; ++s)
                zh_tmp[k * grid_ratio_k + s] = gd.zh[kk] + s * dz;
        }

        // Reconstruct new full level heights.
        z[0] = TF(0.5) * (zh_tmp[0] + zh_tmp[1]);
        for (int k=1; k<ktot_n; ++k)
            z[k] = TF(2) * zh_tmp[k] - z[k-1];

        // Remove top layer from `zh`: it should not be included in the 3D `w` files.
        zh = std::vector<TF>(zh_tmp.begin(), zh_tmp.begin() + ktot_n);
    }
    else
    {
        z  = std::vector<TF>(gd.z.begin()  + gd.kstart, gd.z.begin()  + gd.kend);
        zh = std::vector<TF>(gd.zh.begin() + gd.kstart, gd.zh.begin() + gd.kend);
    }


    auto setup_edge = [&](
        const std::vector<TF>& x_lbc,
        const std::vector<TF>& y_lbc,
        const std::vector<TF>& z_lbc,
        const std::vector<TF>& x,
        const std::vector<TF>& y,
        const std::vector<TF>& z)
    {
        // Save all height levels.
        const int kstart = gd.kgc;
        const int ktot = gd.ktot;

        return NN_interpolator<TF>(x_lbc, y_lbc, z_lbc, x, y, z, gd, md);
    };

    // Scalars.
    for (auto& fld : fields.sp)
    {
        lbc_w.emplace(fld.first, setup_edge(x_w,  y_ew, z, gd.x, gd.y, gd.z));
        lbc_e.emplace(fld.first, setup_edge(x_e,  y_ew, z, gd.x, gd.y, gd.z));
        lbc_s.emplace(fld.first, setup_edge(x_ns, y_s,  z, gd.x, gd.y, gd.z));
        lbc_n.emplace(fld.first, setup_edge(x_ns, y_n,  z, gd.x, gd.y, gd.z));
    }

    // Velocity components.
    lbc_w.emplace("u", setup_edge(xh_w,  y_ew, z, gd.xh, gd.y, gd.z));
    lbc_e.emplace("u", setup_edge(xh_e,  y_ew, z, gd.xh, gd.y, gd.z));
    lbc_s.emplace("u", setup_edge(xh_ns, y_s,  z, gd.xh, gd.y, gd.z));
    lbc_n.emplace("u", setup_edge(xh_ns, y_n,  z, gd.xh, gd.y, gd.z));

    lbc_w.emplace("v", setup_edge(x_w,  yh_ew, z, gd.x, gd.yh, gd.z));
    lbc_e.emplace("v", setup_edge(x_e,  yh_ew, z, gd.x, gd.yh, gd.z));
    lbc_s.emplace("v", setup_edge(x_ns, yh_s,  z, gd.x, gd.yh, gd.z));
    lbc_n.emplace("v", setup_edge(x_ns, yh_n,  z, gd.x, gd.yh, gd.z));

    lbc_w.emplace("w", setup_edge(x_w,  y_ew, zh, gd.x, gd.y, gd.zh));
    lbc_e.emplace("w", setup_edge(x_e,  y_ew, zh, gd.x, gd.y, gd.zh));
    lbc_s.emplace("w", setup_edge(x_ns, y_s,  zh, gd.x, gd.y, gd.zh));
    lbc_n.emplace("w", setup_edge(x_ns, y_n,  zh, gd.x, gd.y, gd.zh));

    if (sw_save_wtop)
    {
        // NN-interpolated `w` at domain top.
        std::vector<TF> x_bc = nnk::arange<TF>(xstart + 0.5 * dx, xend, dx);
        std::vector<TF> y_bc = nnk::arange<TF>(ystart + 0.5 * dy, yend, dy);
        std::vector<TF> zh_bc = {gd.zh[gd.kend]};

        bc_wtop = NN_interpolator<TF>(x_bc, y_bc, zh_bc, gd.x, gd.y, gd.zh, gd, md);
    }

    if (sw_save_buffer)
    {
        // Find start indexes of buffer layer in new vertical coordinates.
        auto it  = std::lower_bound(z.begin(),  z.end(),  zstart_buffer);
        auto ith = std::lower_bound(zh.begin(), zh.end(), zstart_buffer);

        if (it == z.end())
            throw std::runtime_error("Output buffer subdomain is too close to the model top.");

        const int kstart  = std::distance(z.begin(), it);
        const int kstarth = std::distance(zh.begin(), ith);

        std::vector<TF> x_bc = nnk::arange<TF>(xstart + 0.5 * dx, xend, dx);
        std::vector<TF> y_bc = nnk::arange<TF>(ystart + 0.5 * dy, yend, dy);

        std::vector<TF> z_bc(z.begin() + kstart, z.end());
        std::vector<TF> zh_bc(zh.begin() + kstarth, zh.end());

        bc_buffer  = NN_interpolator<TF>(x_bc, y_bc, z_bc,  gd.x, gd.y, gd.z, gd, md);
        bc_bufferh = NN_interpolator<TF>(x_bc, y_bc, zh_bc, gd.x, gd.y, gd.zh, gd, md);
    }
}


template<typename TF>
unsigned long Subdomain<TF>::get_time_limit(unsigned long itime)
{
    unsigned long idtlim = Constants::ulhuge;

    if (sw_subdomain)
    {
        const unsigned long ifreq_bcs = convert_to_itime(savetime_bcs);
        idtlim = std::min(idtlim, ifreq_bcs - itime % ifreq_bcs);

        if (sw_save_buffer)
        {
            const unsigned long ifreq_buf = convert_to_itime(savetime_buffer);
            idtlim = std::min(idtlim, ifreq_buf - itime % ifreq_buf);
        }
    }

    return idtlim;
}


template <typename TF>
void Subdomain<TF>::save_bcs(
        Timeloop<TF>& timeloop)
{
    if (!sw_subdomain || timeloop.in_substep() )
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    auto save_binary = [&](
            NN_interpolator<TF>& lbc,
            const std::string& filename)
    {
        #ifdef USEMPI
        /*
        Save binary using MPI I/O hyperslabs.
        */
        MPI_File fh;
        if (MPI_File_open(md.commxy, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
            return 1;

        int err = 0;
        MPI_Offset fileoff = 0;
        char name[] = "native";
        MPI_Datatype subarray = MPI_DATATYPE_NULL;

        if (lbc.has_data)
        {
            int tot_size[3] = {lbc.ktot_g, lbc.jtot_g, lbc.itot_g};
            int sub_size[3] = {lbc.ktot_s, lbc.jtot_s, lbc.itot_s};
            int sub_start[3] = {0, lbc.j_range.first, lbc.i_range.first};
            int count = lbc.ktot_s * lbc.jtot_s * lbc.itot_s;

            //if (filename.find("buffer") != std::string::npos)
            //std::cout << filename << " | mpi=(" << md.mpicoordx << "," << md.mpicoordy << "), tot=(" << tot_size[0] << "," << tot_size[1] << "," << tot_size[2] << "), sub=" << sub_size[0] << "," << sub_size[1] << "," << sub_size[2] << "), start=" << sub_start[0] << "," << sub_start[1] << "," << sub_start[2] << std::endl;

            MPI_Type_create_subarray(3, tot_size, sub_size, sub_start, MPI_ORDER_C, mpi_fp_type<TF>(), &subarray);
            MPI_Type_commit(&subarray);

            err += MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subarray, name, MPI_INFO_NULL);
            err += MPI_File_write_all(fh, lbc.fld.data(), count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
        }
        else
        {
            err += MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), mpi_fp_type<TF>(), name, MPI_INFO_NULL);
            err += MPI_File_write_all(fh, nullptr, 0, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
        }

        err += MPI_File_close(&fh);

        if (lbc.has_data)
            MPI_Type_free(&subarray);

        return err;

        #else
        /*
        Save binary using simple serial I/O.
        */
        FILE *pFile;
        pFile = fopen(filename.c_str(), "wb");

        if (pFile == NULL)
            return 1;

        fwrite(lbc.fld.data(), sizeof(TF), lbc.fld.size(), pFile);
        fclose(pFile);

        #endif

        return 0;
    };


    if (timeloop.get_itime() % convert_to_itime(savetime_bcs) == 0)
    {
        // Save boundary conditions (lateral + top).
        std::string msg = "Saving subdomain LBCs for time " + std::to_string(timeloop.get_time()) + " ...";
        master.print_message(msg);

        auto process_lbc = [&](
            NN_interpolator<TF>& lbc,
            std::vector<TF>& fld,
            const std::string& name,
            const std::string& loc)
        {
            nnk::nn_interpolate(
                lbc.fld.data(),
                fld.data(),
                lbc.nn_i.data(),
                lbc.nn_j.data(),
                lbc.nn_k.data(),
                lbc.itot_s,
                lbc.jtot_s,
                lbc.ktot_s,
                lbc.istride,
                lbc.jstride,
                lbc.kstride,
                gd.istride,
                gd.jstride,
                gd.kstride);

            // Setup filename with time.
            std::string base_name = "lbc_" + name + "_" + loc + "_out";
            std::string file_name = timeloop.get_io_filename(base_name);

            const int err = save_binary(lbc, file_name);

            if (err > 0)
                throw std::runtime_error("Error saving LBCs.");
        };

        // Save all prognostic fields.
        for (auto& fld : fields.ap)
        {
            process_lbc(lbc_w.at(fld.first), fld.second->fld, fld.first, "west");
            process_lbc(lbc_e.at(fld.first), fld.second->fld, fld.first, "east");
            process_lbc(lbc_s.at(fld.first), fld.second->fld, fld.first, "south");
            process_lbc(lbc_n.at(fld.first), fld.second->fld, fld.first, "north");
        }

        if (sw_save_wtop)
        {
            const int kstart = gd.kend;
            const int ktot = 1;

            nnk::nn_interpolate(
                bc_wtop.fld.data(),
                fields.ap.at("w")->fld.data(),
                bc_wtop.nn_i.data(),
                bc_wtop.nn_j.data(),
                bc_wtop.nn_k.data(),
                bc_wtop.itot_s,
                bc_wtop.jtot_s,
                bc_wtop.ktot_s,
                bc_wtop.istride,
                bc_wtop.jstride,
                bc_wtop.kstride,
                gd.istride,
                gd.jstride,
                gd.kstride);

            // Setup filename with time.
            std::string base_name = "w_top_out";
            std::string file_name = timeloop.get_io_filename(base_name);

            const int err = save_binary(bc_wtop, file_name);

            if (err > 0)
                throw std::runtime_error("Error saving LBCs.");
        }

    }

    if (sw_save_buffer && timeloop.get_itime() % convert_to_itime(savetime_buffer) == 0)
    {
        // Save buffer layer.
        std::string msg = "Saving sub-domain buffer for time " + std::to_string(timeloop.get_time()) + " ...";
        master.print_message(msg);

        for (auto& fld : fields.ap)
        {
            // Switch between full and half levels.
            NN_interpolator<TF>& bc_buff = fld.first == "w" ? bc_bufferh : bc_buffer;

            nnk::nn_interpolate(
                bc_buff.fld.data(),
                fld.second->fld.data(),
                bc_buff.nn_i.data(),
                bc_buff.nn_j.data(),
                bc_buff.nn_k.data(),
                bc_buff.itot_s,
                bc_buff.jtot_s,
                bc_buff.ktot_s,
                bc_buff.istride,
                bc_buff.jstride,
                bc_buff.kstride,
                gd.istride,
                gd.jstride,
                gd.kstride);

            // Setup filename with time.
            std::string base_name = fld.first + "_buffer_out";
            std::string file_name = timeloop.get_io_filename(base_name);

            const int err = save_binary(bc_buff, file_name);

            if (err > 0)
                throw std::runtime_error("Error saving LBCs.");
        }
    }
}

#ifdef FLOAT_SINGLE
template class Subdomain<float>;
#else
template class Subdomain<double>;
#endif
