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

#ifndef FIELD3D
#define FIELD3D

#include <string>

class Master;
class Grid;

class Field3d
{
    public:
        // functions
        Field3d(Grid*, Master*, std::string, std::string, std::string);
        ~Field3d();

        int init();
        // int checkfornan();

        // variables at CPU
        double* data;
        double* databot;
        double* datatop;
        double* datamean;
        double* datagradbot;
        double* datagradtop;
        double* datafluxbot;
        double* datafluxtop;
        std::string name;
        std::string unit;
        std::string longname;
        double visc;

        // Device functions and variables
        void init_device();  ///< Allocate Field3D fields at device 
        void clear_device(); ///< Deallocate Field3D fields at device 

        double* data_g;
        double* databot_g;
        double* datatop_g;
        double* datamean_g;
        double* datagradbot_g;
        double* datagradtop_g;
        double* datafluxbot_g;
        double* datafluxtop_g;

    private:
        Grid* grid;
        Master* master;
};
#endif

