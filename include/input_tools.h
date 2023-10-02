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

#ifndef INPUT_TOOLS_H
#define INPUT_TOOLS_H

#include "master.h"

namespace Input_tools
{
    template<typename T>
    inline void check_item(const T& t) {}

    template<>
    inline void check_item(const std::string& s)
    {
        // Check whether string is empty.
        if (s.empty())
            throw std::runtime_error("Illegal string");
    }

    template<typename T>
    inline T get_item_from_stream(std::istringstream& ss)
    {
        // Read the item from the stringstream, operator >> trims automatically.
        T item;
        if (!(ss >> item))
            throw std::runtime_error("Item does not match type");

        // Check whether stringstream is empty, if not type is incorrect.
        std::string dummy;
        if (ss >> dummy)
            throw std::runtime_error("Item does not match type");

        return item;
    }

    // In case of string, just return the content of the stringstream.
    template<>
    inline std::string get_item_from_stream(std::istringstream& ss)
    {
        std::string item_string = ss.str();
        boost::trim(item_string);
        return item_string;
    }

    // In case of bool, capture string labels true and false.
    template<>
    inline bool get_item_from_stream(std::istringstream& ss)
    {
        bool item;

        std::string item_string = ss.str();
        boost::trim(item_string);
        boost::to_lower(item_string);

        if ( (item_string == "1") || (item_string == "true") )
            item = true;
        else if ( (item_string == "0") || (item_string == "false") )
            item = false;
        else
            throw std::runtime_error("Item does not match type");

        return item;
    }

    inline bool get_line_from_input(std::ifstream& infile, std::string& line, Master& master)
    {
        int has_line = false;
        if (master.get_mpiid() == 0)
        {
            if (std::getline(infile, line))
                has_line = true;
        }

        master.broadcast(&has_line, 1);
        if (has_line)
        {
            // Broadcasting a std::string. This is ugly!
            int line_size = line.size();
            master.broadcast(&line_size, 1);
            if (master.get_mpiid() != 0)
                line.resize(line_size);
            master.broadcast(const_cast<char*>(line.data()), line_size);
        }
        return has_line;
    }
}
#endif
