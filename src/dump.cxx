/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#include <cstdio>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "dump.h"
#include "timeloop.h"
#include "constants.h"
#include "defines.h"


namespace
{
    template<typename TF>
    void remove_ghost_cells(
        TF* const restrict fld_out,
        const TF* const restrict fld_in,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int kstart, const int kend,
        const int igc, const int jgc, const int kgc,
        const int icells, const int jcells)
        {
            const int jstride_in = icells;
            const int kstride_in = icells * jcells;

            const int jstride_out = (icells - 2*igc);
            const int kstride_out = (icells - 2*igc) * (jcells - 2*jgc);

            for (int k=kstart; k<kend; k++)
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk_in = i + j*jstride_in + k*kstride_in;
                        const int ijk_out = (i-igc) + (j-jgc)*jstride_out + (k-kgc)*kstride_out;

                        fld_out[ijk_out] = fld_in[ijk_in];
                    }
        }
}


template<typename TF>
Dump<TF>::Dump(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin):
    master(masterin), grid(gridin), fields(fieldsin),
    field3d_io(master, grid)
{
    swdump = inputin.get_item<bool>("dump", "swdump", "", false);
    swdump_sub = inputin.get_item<bool>("dump", "swdump_sub", "", false);

    if (swdump || swdump_sub)
    {
        // Get the time at which the dump sections are triggered.
        sampletime = inputin.get_item<double>("dump", "sampletime", "");

        // Get list of dump variables.
        dumplist = inputin.get_list<std::string>("dump", "dumplist", "", std::vector<std::string>());

        // Whether to do two consecutive dumps in time
        swdoubledump = inputin.get_item<bool>("dump", "swdoubledump", "", false);
        if (swdoubledump && sampletime != inputin.get_item<double>("time", "savetime", ""))
        {
            std::string msg = "Double dump only works if sampletime is equal to savetime";
            throw std::runtime_error(msg);
        }

        if (swdump_sub)
        {
            mpicoordx = inputin.get_list<int>("dump", "mpicoordx", "", std::vector<int>());
            mpicoordy = inputin.get_list<int>("dump", "mpicoordy", "", std::vector<int>());
        }

        // Crash on empty list.
        if (dumplist.empty())
        {
            std::string msg = "Empty Dump list";
            throw std::runtime_error(msg);
        }
    }
    else
    {
        inputin.flag_as_used("dump", "dumplist", "");
        inputin.flag_as_used("dump", "sampletime", "");
    }

}

template<typename TF>
Dump<TF>::~Dump()
{
}

template<typename TF>
void Dump<TF>::init()
{
    if (!(swdump || swdump_sub))
        return;

    isampletime = convert_to_itime(sampletime);
}

template<typename TF>
void Dump<TF>::create()
{
    /* All classes (fields, thermo) have removed their dump-variables from
       dumplist by now. If it isn't empty, print warnings for invalid variables */
    if (!dumplist.empty())
    {
        for (auto& it : dumplist)
            master.print_warning("field %s in [dump][dumplist] is illegal\n", it.c_str());
    }
}

template<typename TF>
unsigned long Dump<TF>::get_time_limit(unsigned long itime)
{
    if (!(swdump || swdump_sub))
        return Constants::ulhuge;

    return isampletime - itime % isampletime;
}

template<typename TF>
bool Dump<TF>::do_dump(unsigned long itime, unsigned long idt)
{
    // Check if dump is enabled.
    if (!(swdump || swdump_sub))
        return false;

    // Check if current time step is dump time.
    if (itime % isampletime != 0)
    {
        if (((itime + idt) % isampletime == 0) && swdoubledump)
            return true;
        else
            return false;
    }

    // Return true such that dumps are created
    return true;
}

template<typename TF>
std::vector<std::string>& Dump<TF>::get_dumplist()
{
    return dumplist;
}

template<typename TF>
void Dump<TF>::save_dump(TF* data, const std::string& varname, int iotime)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const double no_offset = 0.;
    char filename[256];

    if (swdump)
    {
        std::snprintf(filename, 256, "%s.%07d", varname.c_str(), iotime);
        std::ifstream infile(filename);

        if (infile.good())
        {
            master.print_message("%s already exists\n", filename);
        }
        else
        {
            auto tmp1 = fields.get_tmp();
            auto tmp2 = fields.get_tmp();

            if (field3d_io.save_field3d(
                        data,
                        tmp1->fld.data(),
                        tmp2->fld.data(),
                        filename,
                        no_offset,
                        gd.kstart,
                        gd.kend))
            {
                master.print_message("Saving \"%s\" ... FAILED\n", filename);
                throw std::runtime_error("Writing error in dump");
            }

            fields.release_tmp(tmp1);
            fields.release_tmp(tmp2);
        }
    }

    if (swdump_sub)
    {
        // Check combinations. This way you don't need to specify each combination manually.
        bool save_dump = false;
        auto it_x = std::find(mpicoordx.begin(), mpicoordx.end(), md.mpicoordx);
        auto it_y = std::find(mpicoordy.begin(), mpicoordy.end(), md.mpicoordy);

        if (it_x != mpicoordx.end() && it_y != mpicoordy.end())
            save_dump = true;

        if (save_dump)
        {
            auto data_nogc = fields.get_tmp();

            remove_ghost_cells(
                data_nogc->fld.data(),
                data,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.igc, gd.jgc, gd.kgc,
                gd.icells, gd.jcells);

            std::snprintf(filename, 256, "%s.%04d.%04d.%07d", varname.c_str(), md.mpicoordx, md.mpicoordy, iotime);
            std::ifstream infile(filename);

            FILE *pFile;
            pFile = fopen(filename, "wb");

            if (pFile == NULL)
            {
                #ifdef USEMPI
                std::cout << "SINGLE PROCESS EXCEPTION: saving dump failed." << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
                #else
                throw std::runtime_error("Saving dump failed.");
                #endif
            }

            fwrite(data_nogc->fld.data(), sizeof(TF), gd.imax*gd.jmax*gd.kmax, pFile);
            fclose(pFile);

            fields.release_tmp(data_nogc);
        }
    }
}


#ifdef FLOAT_SINGLE
template class Dump<float>;
#else
template class Dump<double>;
#endif
