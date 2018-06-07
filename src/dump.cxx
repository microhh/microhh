/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "dump.h"
//#include "thermo.h"
#include "timeloop.h"
#include "constants.h"
#include "defines.h"

template<typename TF>
Dump<TF>::Dump(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin):
    master(masterin), grid(gridin), fields(fieldsin),
    field3d_io(master, grid)
{
    swdump = inputin.get_item<bool>("dump", "swdump", "", false);

    if (swdump)
    {
       // Get the time at which the dump sections are triggered.
        sampletime = inputin.get_item<double>("dump", "sampletime", "");

        // Get list of dump variables.
        dumplist = inputin.get_list<std::string>("dump", "dumplist", "", std::vector<std::string>());

        // Crash on empty list.
        if (dumplist.empty())
        {
            master.print_error("empty dump list\n");
            throw 1;
        }
    }
}

template<typename TF>
Dump<TF>::~Dump()
{
}

template<typename TF>
void Dump<TF>::init(double ifactor)
{
    if (!swdump)
        return;

    isampletime = static_cast<unsigned long>(ifactor * sampletime);

    field3d_io.init();
}

template<typename TF>
void Dump<TF>::create()
{
    /* All classes (fields, thermo) have removed their dump-variables from
       dumplist by now. If it isnt empty, print warnings for invalid variables */
    if (!dumplist.empty())
    {
        for (auto& it : dumplist)
            master.print_warning("field %s in [dump][dumplist] is illegal\n", it.c_str());
    }
}

template<typename TF>
unsigned long Dump<TF>::get_time_limit(unsigned long itime)
{
    if (!swdump)
        return Constants::ulhuge;

    return isampletime - itime % isampletime;
}

template<typename TF>
bool Dump<TF>::do_dump(unsigned long itime)
{
    // check if dump are enabled
    if (!swdump)
        return false;
    // check if time for execution
    if (itime % isampletime != 0)
        return false;
    // return true such that column are computed
    return true;
}

template<typename TF>
std::vector<std::string>* Dump<TF>::get_dumplist()
{
    return &dumplist;
}

template<typename TF>
void Dump<TF>::save_dump(TF* data, std::string varname, int iotime)
{
    const double no_offset = 0.;
    char filename[256];

    std::sprintf(filename, "%s.%07d", varname.c_str(), iotime);
    std::ifstream infile(filename);
    if(infile.good())
    {
        master.print_message("%s already exists\n", filename);
    }
    else
    {
        master.print_message("Saving \"%s\" ... ", filename);

        auto tmpfld1 = fields.get_tmp();
        auto tmpfld2 = fields.get_tmp();
        auto tmp1 = tmpfld1.get();
        auto tmp2 = tmpfld1.get();

        if (field3d_io.save_field3d(data, tmp1->fld.data(), tmp2->fld.data(), filename, no_offset))
        {
            master.print_message("FAILED\n");
            throw 1;
        }
        else
        {
            master.print_message("OK\n");
        }
        fields.release_tmp(tmpfld1);
        fields.release_tmp(tmpfld2);
    }
}


template class Dump<double>;
template class Dump<float>;
