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

#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <sstream>
#include "master.h"

void Master::print_message(const char *format, ...)
{
    if (md.mpiid == 0)
    {
        va_list args;
        va_start(args, format);
        std::vfprintf(stdout, format, args);
        va_end(args);
    }
}

void Master::print_message(const std::ostringstream& ss)
{
    if (md.mpiid == 0)
        std::cout << ss.str();
}

void Master::print_message(const std::string& s)
{
    if (md.mpiid == 0)
        std::cout << s << std::endl;
}

void Master::print_warning(const char *format, ...)
{
    std::string warningstr("WARNING: ");
    warningstr += std::string(format);

    const char *warningformat = warningstr.c_str();

    if (md.mpiid == 0)
    {
        va_list args;
        va_start(args, format);
        std::vfprintf(stdout, warningformat, args);
        va_end(args);
    }
}

void Master::print_warning(const std::ostringstream& ss)
{
    if (md.mpiid == 0)
        std::cout << "WARNING: " << ss.str();
}

void Master::print_warning(const std::string& s)
{
    if (md.mpiid == 0)
        std::cout << "WARNING: " << s << std::endl;
}

bool Master::at_wall_clock_limit()
{
    const double wall_clock_time_left = wall_clock_end - get_wall_clock_time();
    const double ten_minutes = 10.*60.;

    if (wall_clock_time_left < ten_minutes)
        return true;
    else
        return false;
}
