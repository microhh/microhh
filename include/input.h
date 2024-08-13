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

#ifndef INPUT_H
#define INPUT_H

#include <map>
#include <vector>
#include <string>

class Master;

class Input
{
    public:
        Input(Master&, const std::string&);
        template<typename T> T get_item(const std::string&, const std::string&, const std::string&);
        template<typename T> T get_item(const std::string&, const std::string&, const std::string&, const T);
        template<typename T> std::vector<T> get_list(const std::string&, const std::string&, const std::string&);
        template<typename T> std::vector<T> get_list(const std::string&, const std::string&, const std::string&, const std::vector<T>);
        // void print_itemlist();
        void print_unused_items();
        void flag_as_used(const std::string&, const std::string&, const std::string&);

        typedef std::map<std::string, std::map< std::string, std::map<std::string, std::pair<std::string, bool>>>> Itemlist;

    private:
        Master& master;
        Itemlist itemlist;
};
#endif
