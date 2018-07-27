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


#include <netcdf>
using namespace netCDF;

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;


template<typename TF>
class Column
{
    public:
        Column(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Column();

        void init(double);
        void create(Input&, int, std::string);

        unsigned long get_time_limit(unsigned long);
        bool get_switch() { return swcolumn; }
        void exec(int, double, unsigned long);
        bool do_column(unsigned long);


        // Interface functions.
        void add_prof(std::string, std::string, std::string, std::string);
        void calc_column(std::string, const TF* const,
                       const TF);


    private:

        // struct for profiles
        struct Prof_var
        {
            NcVar ncvar;
            std::vector<TF> data;
        };

        using Prof_map = std::map<std::string, Prof_var>;

        // structure
        struct Column_struct
        {
            std::vector<int> coord;
            NcFile* data_file;
            NcDim z_dim;
            NcDim zh_dim;
            NcDim t_dim;
            NcVar iter_var;
            NcVar t_var;
            Prof_map profs;
        };

        std::vector<Column_struct> columns;

    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool swcolumn;           ///< Statistics on/off switch

        int statistics_counter;
        double sampletime;
        unsigned long isampletime;

};
#endif
