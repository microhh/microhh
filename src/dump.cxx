/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "dump.h"
#include "model.h"
#include "thermo.h"
#include "timeloop.h"
#include "constants.h"
#include "defines.h"

Dump::Dump(Model *modelin, Input *inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    int nerror = 0;
    nerror += inputin->get_item(&swdump, "dump", "swdump", "", "0");

    if (swdump == "1")
    {  
        nerror += inputin->get_item(&sampletime, "dump", "sampletime", "");
        nerror += inputin->get_list(&dumplist ,  "dump", "dumplist" ,  "");
    }  

    if (nerror)
        throw 1;
}

Dump::~Dump()
{
}

void Dump::init(double ifactor)
{
    if (swdump == "0")
        return;

    isampletime = (unsigned long)(ifactor * sampletime);
}

void Dump::create()
{  
    /* All classes (fields, thermo) have removed their dump-variables from
       dumplist by now. If it isnt empty, print warnings for invalid variables */
    if (dumplist.size() > 0)
    {
        for (std::vector<std::string>::const_iterator it=dumplist.begin(); it!=dumplist.end(); ++it)
            master->print_warning("field %s in [dump][dumplist] is illegal\n", it->c_str());
    } 
}

unsigned long Dump::get_time_limit(unsigned long itime)
{
    if (swdump == "0")
        return Constants::ulhuge;

    return isampletime - itime % isampletime;
}

std::string Dump::get_switch()
{
    return swdump;
}

std::vector<std::string>* Dump::get_dumplist()
{
    return &dumplist;
}

bool Dump::do_dump()
{
    if (swdump == "0")
        return false;

    if (model->timeloop->get_itime() % isampletime == 0)
        return true;
    else
        return false;
}

void Dump::save_dump(double * restrict data, double * restrict tmp, std::string varname)
{
    const double NoOffset = 0.;
    char filename[256];

    std::sprintf(filename, "%s.%07d", varname.c_str(), model->timeloop->get_iotime());
    master->print_message("Saving \"%s\" ... ", filename);

    if (grid->save_field3d(data, tmp, fields->atmp["tmp2"]->data, filename, NoOffset))
    {
        master->print_message("FAILED\n");
        throw 1;
    }  
    else
    {
        master->print_message("OK\n");
    }
}
