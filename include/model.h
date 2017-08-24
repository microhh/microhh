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

#ifndef MODEL
#define MODEL

#include <string>

class Master;
class Input;
class Grid;

class Model
{
    public:
        Model(Master*, Input*);
        ~Model();

        void init();
        void load();
        void save();
        void exec();

    private:
        Master* master;
        Input*  input;
        Grid*   grid;

        void delete_objects();
        void print_status();
        void calc_stats(std::string);
        void set_time_step();
};
#endif
