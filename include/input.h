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

#ifndef INPUT
#define INPUT

#include <map>
#include <string>
#include <vector>

// Forward declaration to avoid circular dependency.
class Master;

typedef std::map<std::string, std::vector<double> > Data_map;

class Input
{
    public:
        Input(Master*);
        ~Input();

        void clear();

        // Item retrieval functions
        int get_item(int*        , std::string, std::string, std::string);
        int get_item(int*        , std::string, std::string, std::string, int);
        int get_item(double*     , std::string, std::string, std::string);
        int get_item(double*     , std::string, std::string, std::string, double);
        int get_item(bool*       , std::string, std::string, std::string);
        int get_item(bool*       , std::string, std::string, std::string, bool);
        int get_item(std::string*, std::string, std::string, std::string);
        int get_item(std::string*, std::string, std::string, std::string, std::string);

        // List retrieval functions
        int get_list(std::vector<int> *        , std::string, std::string, std::string);
        int get_list(std::vector<double> *     , std::string, std::string, std::string);
        int get_list(std::vector<std::string> *, std::string, std::string, std::string);

        int get_prof(double*, std::string, int size);
        int get_time(double**, std::vector<double>*, std::string);
        int get_time_prof(double**, std::vector<double>*, std::string, int);

        void print_unused();
        void flag_as_used(std::string, std::string);

    private:
        Master* master;

        int read_ini_file();
        int read_data_file(Data_map*, std::string, bool);

        template <class valuetype>
        int parse_item(valuetype*, std::string, std::string, std::string, bool, valuetype);

        template <class valuetype>
        int parse_list(std::vector<valuetype> *, std::string, std::string, std::string);

        int check_item_exists(std::string, std::string, std::string el="default");
        int check_item(int*        , std::string, std::string, std::string el="default");
        int check_item(double*     , std::string, std::string, std::string el="default");
        int check_item(bool*       , std::string, std::string, std::string el="default");
        int check_item(std::string*, std::string, std::string, std::string el="default");

        // list retrieval
        int check_list(std::vector<int>*        , std::string, std::string, std::string el="default");
        int check_list(std::vector<double>*     , std::string, std::string, std::string el="default");
        int check_list(std::vector<std::string>*, std::string, std::string, std::string el="default");

        struct Input_type
        {
            std::string data;
            bool isused;
        };
        typedef std::map<std::string, Input_type  > Input_map_1d;
        typedef std::map<std::string, Input_map_1d> Input_map_2d;
        typedef std::map<std::string, Input_map_2d> Input_map;

        Input_map inputlist;
        Data_map proflist;
        Data_map timelist;

        std::string isused;
};
#endif
