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

#ifndef COLUMN
#define COLUMN

//#include <netcdfcpp.h>
#include <netcdf>
using namespace netCDF;

class Master;
class Model;
class Grid;
class Fields;

// struct for profiles
struct Column_var
{
    NcVar ncvar;
    double* data;
};

// typedefs for containers of profiles and time series
typedef std::map<std::string, Column_var> Column_map;

class Column
{
    public:
        Column(Model*, Input*);
        ~Column();

        void init(double);
        void create(int);

        unsigned long get_time_limit(unsigned long);
        void exec(int, double, unsigned long);
        bool doColumn();
        std::string get_switch();

        // Interface functions.
        void add_prof(std::string, std::string, std::string, std::string);
        void calc_column(double* const, const double* const,
                       const double);
        std::string name;
        NcFile* dataFile;
        NcDim z_dim;
        NcDim zh_dim;
        NcDim t_dim;
        NcVar iter_var;
        NcVar t_var;
        Column_map profs;

    private:
        int ncolumn;

        // mask calculations
        void calc_column(double* const, const double* const,
                       const double, const int[2]);

    protected:
        Model*  model;
        Grid*   grid;
        Fields* fields;
        Master* master;

        double sampletime;
        long isampletime;

        std::string swcolumn;

};
#endif
