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

#ifndef FIELD3D
#define FIELD3D

#include <string>

class Master;
template<typename> class Grid;

template<typename TF>
class Field3d
{
    public:
        // functions
        Field3d(Master&, Grid<TF>&, std::string, std::string, std::string);
        ~Field3d();

        int init();

        // variables at CPU
        std::vector<TF> data;
        std::vector<TF> databot;
        std::vector<TF> datatop;
        std::vector<TF> datamean;
        std::vector<TF> datagradbot;
        std::vector<TF> datagradtop;
        std::vector<TF> datafluxbot;
        std::vector<TF> datafluxtop;

        std::string name;
        std::string unit;
        std::string longname;

        TF visc;

        /*
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
        */

    private:
        Master& master;
        Grid<TF>& grid;
};
#endif

