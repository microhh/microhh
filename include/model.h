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

#ifndef MODEL
#define MODEL

#include <string>

class Master;
class Input;
class Grid;
class Fields;
class Boundary;
class Timeloop;
class Advec;
class Diff;
class Pres;
class Force;
class Thermo;
class Buffer;
class Stats;
class Cross;
class Dump;
class Budget;

class Model
{
    public:
        Model(Master*, Input*);
        ~Model();

        void init();
        void load();
        void save();
        void exec();

        // Make the pointers public for use in other classes.
        // TODO maybe it is safer to create get functions
        Master* master;
        Input*  input;
        Grid*   grid;
        Fields* fields;

        // Model operators.
        Boundary* boundary;
        Timeloop* timeloop;
        Advec*    advec;
        Diff*     diff;
        Pres*     pres;  
        Force*    force;   
        Thermo*   thermo;
        Buffer*   buffer;

        // Postprocessing and output modules.
        Stats*  stats;
        Cross*  cross;
        Dump*   dump;
        Budget* budget;

    private:
        // list of masks for statistics
        std::vector<std::string> masklist;

        void delete_objects();

        void print_status();
        void calc_stats(std::string);
        void set_time_step();
};
#endif
