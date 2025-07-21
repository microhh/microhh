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

#include "subdomain.h"

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "constants.h"

namespace blk = boundary_lateral_kernels;

namespace
{
    #ifdef USEMPI
    template<typename TF> MPI_Datatype mpi_fp_type();
    template<> MPI_Datatype mpi_fp_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_fp_type<float>() { return MPI_FLOAT; }
    #endif
}

template<typename TF>
Subdomain<TF>::Subdomain(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin), field3d_io(masterin, gridin)
{
    sw_subdomain = inputin.get_item<bool>("subdomain", "sw_subdomain", "", false);


    //if (sw_subdomain)
    //{
    //    sw_save_wtop = inputin.get_item<bool>("subdomain", "sw_save_wtop", "", false);
    //    sw_save_buffer = inputin.get_item<bool>("subdomain", "sw_save_buffer", "", false);

    //    xstart_sub = inputin.get_item<TF>("subdomain", "xstart", "");
    //    xend_sub   = inputin.get_item<TF>("subdomain", "xend", "");
    //    ystart_sub = inputin.get_item<TF>("subdomain", "ystart", "");
    //    yend_sub   = inputin.get_item<TF>("subdomain", "yend", "");

    //    grid_ratio_sub = inputin.get_item<int>("subdomain", "grid_ratio", "");
    //    n_ghost_sub = inputin.get_item<int>("subdomain", "n_ghost", "");
    //    n_sponge_sub = inputin.get_item<int>("subdomain", "n_sponge", "");

    //    savetime_bcs = inputin.get_item<int>("subdomain", "savetime_bcs", "");

    //    if (sw_save_buffer)
    //    {
    //        savetime_buffer = inputin.get_item<int>("subdomain", "savetime_buffer", "");
    //        zstart_buffer = inputin.get_item<TF>("subdomain", "zstart_buffer", "");
    //    }
    //}
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

//    bool error = false;
//
//    auto check_subdomain = [&](const TF value, const TF spacing, const std::string& name)
//    {
//        if (!is_equal(std::fmod(value, spacing), TF(0)))
//        {
//            error = true;
//            master.print_message(
//                "ERROR: " + name + " must be an integer multiple of " + std::to_string(spacing) + ".");
//        }
//    };
//
//    // Sub-domain has perfectly align with parent grid at x/y half levels.
//    check_subdomain(xstart_sub, gd.dx, "xstart_sub");
//    check_subdomain(xend_sub, gd.dx, "xend_sub");
//    check_subdomain(ystart_sub, gd.dy, "ystart_sub");
//    check_subdomain(yend_sub, gd.dy, "yend_sub");
//
//    if (error)
//        throw std::runtime_error("Sub-domain boundaries not aligned with parent grid.");
//
//    // Initialise LBCs instance.
//    xsize_sub = xend_sub - xstart_sub;
//    ysize_sub = yend_sub - ystart_sub;
//
//    itot_sub = static_cast<int>(xsize_sub / gd.dx * grid_ratio_sub + 0.5);
//    jtot_sub = static_cast<int>(ysize_sub / gd.dy * grid_ratio_sub + 0.5);
//
//    dx_sub = xsize_sub / itot_sub;
//    dy_sub = ysize_sub / jtot_sub;
//
//    // Coordinates of LBCs at different grid locations.
//    // Full levels x/y.
//    std::vector<TF> x_ns = blk::arange<TF>(xstart_sub - (n_ghost_sub  - 0.5) * dx_sub, xend_sub   + n_ghost_sub  * dx_sub, dx_sub);
//    std::vector<TF> x_w  = blk::arange<TF>(xstart_sub - (n_ghost_sub  - 0.5) * dx_sub, xstart_sub + n_sponge_sub * dx_sub, dx_sub);
//    std::vector<TF> x_e  = blk::arange<TF>(xend_sub   - (n_sponge_sub - 0.5) * dx_sub, xend_sub   + n_ghost_sub  * dx_sub, dx_sub);
//
//    std::vector<TF> y_ew = blk::arange<TF>(ystart_sub - (n_ghost_sub  - 0.5) * dy_sub, yend_sub   + n_ghost_sub  * dy_sub, dy_sub);
//    std::vector<TF> y_n  = blk::arange<TF>(yend_sub   - (n_sponge_sub - 0.5) * dy_sub, yend_sub   + n_ghost_sub  * dy_sub, dy_sub);
//    std::vector<TF> y_s  = blk::arange<TF>(ystart_sub - (n_ghost_sub  - 0.5) * dy_sub, ystart_sub + n_sponge_sub * dy_sub, dy_sub);
//
//    // Half level `x` for `u`.
//    std::vector<TF> xh_ns = blk::arange<TF>(xstart_sub - n_ghost_sub  * dx_sub, xend_sub   + n_ghost_sub  * dx_sub,            dx_sub);
//    std::vector<TF> xh_w  = blk::arange<TF>(xstart_sub - n_ghost_sub  * dx_sub, xstart_sub + n_sponge_sub * dx_sub + dx_sub/2, dx_sub);
//    std::vector<TF> xh_e  = blk::arange<TF>(xend_sub   - n_sponge_sub * dx_sub, xend_sub   + n_ghost_sub  * dx_sub,            dx_sub);
//
//    // Half level `y` for `v`.
//    std::vector<TF> yh_ew = blk::arange<TF>(ystart_sub - n_ghost_sub  * dy_sub, yend_sub   + n_ghost_sub  * dy_sub,            dy_sub);
//    std::vector<TF> yh_n  = blk::arange<TF>(yend_sub   - n_sponge_sub * dy_sub, yend_sub   + n_ghost_sub  * dy_sub,            dy_sub);
//    std::vector<TF> yh_s  = blk::arange<TF>(ystart_sub - n_ghost_sub  * dy_sub, ystart_sub + n_sponge_sub * dy_sub + dy_sub/2, dy_sub);
//
//    auto setup_edge = [&](
//        const std::vector<TF>& x_lbc,
//        const std::vector<TF>& y_lbc,
//        const std::vector<TF>& x,
//        const std::vector<TF>& y)
//    {
//        return Nn_data<TF>(x_lbc, y_lbc, x, y, gd, md, gd.ktot);
//    };
//
//    // Scalars.
//    for (auto& fld : fields.sp)
//    {
//        lbc_sub_w.emplace(fld.first, setup_edge(x_w,  y_ew, gd.x, gd.y));
//        lbc_sub_e.emplace(fld.first, setup_edge(x_e,  y_ew, gd.x, gd.y));
//        lbc_sub_s.emplace(fld.first, setup_edge(x_ns, y_s,  gd.x, gd.y));
//        lbc_sub_n.emplace(fld.first, setup_edge(x_ns, y_n,  gd.x, gd.y));
//    }
//
//    // Velocity components.
//    lbc_sub_w.emplace("u", setup_edge(xh_w,  y_ew, gd.xh, gd.y));
//    lbc_sub_e.emplace("u", setup_edge(xh_e,  y_ew, gd.xh, gd.y));
//    lbc_sub_s.emplace("u", setup_edge(xh_ns, y_s,  gd.xh, gd.y));
//    lbc_sub_n.emplace("u", setup_edge(xh_ns, y_n,  gd.xh, gd.y));
//
//    lbc_sub_w.emplace("v", setup_edge(x_w,  yh_ew, gd.x, gd.yh));
//    lbc_sub_e.emplace("v", setup_edge(x_e,  yh_ew, gd.x, gd.yh));
//    lbc_sub_s.emplace("v", setup_edge(x_ns, yh_s,  gd.x, gd.yh));
//    lbc_sub_n.emplace("v", setup_edge(x_ns, yh_n,  gd.x, gd.yh));
//
//    lbc_sub_w.emplace("w", setup_edge(x_w,  y_ew, gd.x, gd.y));
//    lbc_sub_e.emplace("w", setup_edge(x_e,  y_ew, gd.x, gd.y));
//    lbc_sub_s.emplace("w", setup_edge(x_ns, y_s,  gd.x, gd.y));
//    lbc_sub_n.emplace("w", setup_edge(x_ns, y_n,  gd.x, gd.y));
//
//    if (sw_save_wtop)
//    {
//        // NN-interpolated `w` at domain top.
//        std::vector<TF> x = blk::arange<TF>(xstart_sub + 0.5 * dx_sub, xend_sub, dx_sub);
//        std::vector<TF> y = blk::arange<TF>(ystart_sub + 0.5 * dy_sub, yend_sub, dy_sub);
//
//        bc_wtop_sub = Nn_data<TF>(x, y, gd.x, gd.y, gd, md, 1);
//    }
//
//    if (sw_save_buffer)
//    {
//        // Find start indexes of buffer layer.
//        buffer_kstart  = gd.kstart;
//        buffer_kstarth = gd.kstart;
//
//        for (int k=gd.kstart; k<gd.kend; ++k)
//        {
//            if (gd.z[k] < zstart_buffer)
//                ++buffer_kstart;
//            if (gd.zh[k] < zstart_buffer)
//                ++buffer_kstarth;
//        }
//
//        if (buffer_kstarth == gd.kend)
//            throw std::runtime_error("Output buffer subdomain is too close to the model top.");
//
//        std::vector<TF> x = blk::arange<TF>(xstart_sub + 0.5 * dx_sub, xend_sub, dx_sub);
//        std::vector<TF> y = blk::arange<TF>(ystart_sub + 0.5 * dy_sub, yend_sub, dy_sub);
//
//        const int ksize = gd.kend - buffer_kstart;
//        const int ksizeh = gd.kend - buffer_kstarth;
//
//        bc_buffer  = Nn_data<TF>(x, y, gd.x, gd.y, gd, md, ksize);
//        bc_bufferh = Nn_data<TF>(x, y, gd.x, gd.y, gd, md, ksizeh);
//    }
}


template<typename TF>
unsigned long Subdomain<TF>::get_time_limit(unsigned long itime)
{
    unsigned long idtlim = Constants::ulhuge;

    //if (sw_subdomain)
    //{
    //    const unsigned long ifreq_bcs = convert_to_itime(savetime_bcs);
    //    idtlim = std::min(idtlim, ifreq_bcs - itime % ifreq_bcs);

    //    if (sw_save_buffer)
    //    {
    //        const unsigned long ifreq_buf = convert_to_itime(savetime_buffer);
    //        idtlim = std::min(idtlim, ifreq_buf - itime % ifreq_buf);
    //    }
    //}

    return idtlim;
}


template <typename TF>
void Subdomain<TF>::save_bcs(
        Timeloop<TF>& timeloop)
{
    if (!sw_subdomain || timeloop.in_substep() )
        return;

//    auto& gd = grid.get_grid_data();
//    auto& md = master.get_MPI_data();
//
//    auto save_binary = [&](
//            Nn_data<TF>& lbc,
//            const std::string& filename)
//    {
//        #ifdef USEMPI
//        /*
//        Save binary using MPI I/O hyperslabs.
//        */
//        MPI_File fh;
//        if (MPI_File_open(md.commxy, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
//            return 1;
//
//        int err = 0;
//        MPI_Offset fileoff = 0;
//        char name[] = "native";
//        MPI_Datatype subarray = MPI_DATATYPE_NULL;
//
//        if (lbc.has_data)
//        {
//            int tot_size[3] = {lbc.ktot_g, lbc.jtot_g, lbc.itot_g};
//            int sub_size[3] = {lbc.ktot_s, lbc.jtot_s, lbc.itot_s};
//            int sub_start[3] = {0, lbc.j_range.first, lbc.i_range.first};
//            int count = lbc.ktot_s * lbc.jtot_s * lbc.itot_s;
//
//            //if (filename.find("buffer") != std::string::npos)
//            //    std::cout << filename << " | mpi=(" << md.mpicoordx << "," << md.mpicoordy << "), tot=(" << tot_size[0] << "," << tot_size[1] << "," << tot_size[2] << "), sub=" << sub_size[0] << "," << sub_size[1] << "," << sub_size[2] << "), start=" << sub_start[0] << "," << sub_start[1] << "," << sub_start[2] << std::endl;
//
//            MPI_Type_create_subarray(3, tot_size, sub_size, sub_start, MPI_ORDER_C, mpi_fp_type<TF>(), &subarray);
//            MPI_Type_commit(&subarray);
//
//            err += MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subarray, name, MPI_INFO_NULL);
//            err += MPI_File_write_all(fh, lbc.fld.data(), count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
//        }
//        else
//        {
//            err += MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), mpi_fp_type<TF>(), name, MPI_INFO_NULL);
//            err += MPI_File_write_all(fh, nullptr, 0, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
//        }
//
//        err += MPI_File_close(&fh);
//
//        if (lbc.has_data)
//            MPI_Type_free(&subarray);
//
//        return err;
//
//        #else
//        /*
//        Save binary using simple serial I/O.
//        */
//        FILE *pFile;
//        pFile = fopen(filename.c_str(), "wb");
//
//        if (pFile == NULL)
//            return 1;
//
//        fwrite(lbc.fld.data(), sizeof(TF), lbc.fld.size(), pFile);
//        fclose(pFile);
//
//        #endif
//
//        return 0;
//    };
//
//
//    if (timeloop.get_itime() % convert_to_itime(savetime_bcs) == 0)
//    {
//        // Save boundary conditions (lateral + top).
//        std::string msg = "Saving sub-domain LBCs for time " + std::to_string(timeloop.get_time()) + " ...";
//        master.print_message(msg);
//
//        auto process_lbc = [&](
//            Nn_data<TF>& lbc,
//            std::vector<TF>& fld,
//            const std::string& name,
//            const std::string& loc)
//        {
//            blk::nn_interpolate(
//                lbc.fld.data(),
//                fld.data(),
//                lbc.nn_i.data(),
//                lbc.nn_j.data(),
//                lbc.itot_s,
//                lbc.jtot_s,
//                gd.ktot,
//                gd.kgc,
//                lbc.istride,
//                lbc.jstride,
//                lbc.kstride,
//                gd.istride,
//                gd.jstride,
//                gd.kstride);
//
//            // Setup filename with time.
//            std::string base_name = "lbc_" + name + "_" + loc + "_out";
//            std::string file_name = timeloop.get_io_filename(base_name);
//
//            const int err = save_binary(lbc, file_name);
//
//            if (err > 0)
//                throw std::runtime_error("Error saving LBCs.");
//        };
//
//        // Save all prognostic fields.
//        for (auto& fld : fields.ap)
//        {
//            process_lbc(lbc_sub_w.at(fld.first), fld.second->fld, fld.first, "west");
//            process_lbc(lbc_sub_e.at(fld.first), fld.second->fld, fld.first, "east");
//            process_lbc(lbc_sub_s.at(fld.first), fld.second->fld, fld.first, "south");
//            process_lbc(lbc_sub_n.at(fld.first), fld.second->fld, fld.first, "north");
//        }
//
//        if (sw_save_wtop)
//        {
//            const int kstart = gd.kend;
//            const int ktot = 1;
//
//            blk::nn_interpolate(
//                bc_wtop_sub.fld.data(),
//                fields.ap.at("w")->fld.data(),
//                bc_wtop_sub.nn_i.data(),
//                bc_wtop_sub.nn_j.data(),
//                bc_wtop_sub.itot_s,
//                bc_wtop_sub.jtot_s,
//                ktot,
//                kstart,
//                bc_wtop_sub.istride,
//                bc_wtop_sub.jstride,
//                bc_wtop_sub.kstride,
//                gd.istride,
//                gd.jstride,
//                gd.kstride);
//
//            // Setup filename with time.
//            std::string base_name = "w_top_out";
//            std::string file_name = timeloop.get_io_filename(base_name);
//
//            const int err = save_binary(bc_wtop_sub, file_name);
//
//            if (err > 0)
//                throw std::runtime_error("Error saving LBCs.");
//        }
//
//    }
//
//
//    if (sw_save_buffer && timeloop.get_itime() % convert_to_itime(savetime_buffer) == 0)
//    {
//        // Save buffer layer.
//        std::string msg = "Saving sub-domain buffer for time " + std::to_string(timeloop.get_time()) + " ...";
//        master.print_message(msg);
//
//        for (auto& fld : fields.ap)
//        {
//            const int kstart = fld.first == "w" ? buffer_kstarth : buffer_kstart;
//            Nn_data<TF>& bc_buff = fld.first == "w" ? bc_bufferh : bc_buffer;
//
//            blk::nn_interpolate(
//                bc_buff.fld.data(),
//                fld.second->fld.data(),
//                bc_buff.nn_i.data(),
//                bc_buff.nn_j.data(),
//                bc_buff.itot_s,
//                bc_buff.jtot_s,
//                bc_buff.ktot_s,
//                kstart,
//                bc_buff.istride,
//                bc_buff.jstride,
//                bc_buff.kstride,
//                gd.istride,
//                gd.jstride,
//                gd.kstride);
//
//            // Setup filename with time.
//            std::string base_name = fld.first + "_buffer_out";
//            std::string file_name = timeloop.get_io_filename(base_name);
//
//            const int err = save_binary(bc_buff, file_name);
//
//            if (err > 0)
//                throw std::runtime_error("Error saving LBCs.");
//        }
//    }
}

#ifdef FLOAT_SINGLE
template class Subdomain<float>;
#else
template class Subdomain<double>;
#endif
