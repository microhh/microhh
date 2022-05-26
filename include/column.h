/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef COLUMN_H
#define COLUMN_H

class Master;
class Input;
class Netcdf_file;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;
template<typename> class Netcdf_variable;

template<typename TF>
class Column
{
    public:
        Column(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Column();

        void init(double);
        void create(Input&, Timeloop<TF>&, std::string);

        unsigned long get_time_limit(unsigned long);
        bool get_switch() { return swcolumn; }
        int get_n_columns() {return columns.size(); }
        void get_column_locations(std::vector<int>&, std::vector<int>&);
        void exec(int, double, unsigned long);
        bool do_column(unsigned long);

        // Interface functions.
        void add_prof(std::string, std::string, std::string, std::string);
        void add_time_series(std::string, std::string, std::string);

        void calc_column(std::string, const TF* const, const TF, const bool copy_from_gpu=true);
        void calc_time_series(std::string, const TF* const, const TF);
        void set_individual_column(std::string, const TF*, const TF, const int, const int);

        #ifdef USECUDA
        void prepare_device();
        void clear_device();
        int* get_column_location_g(const std::string&);
        #endif

    private:
        // Struct for profiles.
        struct Prof_var
        {
            Netcdf_variable<TF> ncvar;
            std::vector<TF> data;
        };

        // Struct for time series.
        struct Time_var
        {
            Netcdf_variable<TF> ncvar;
            TF data;
        };

        using Prof_map = std::map<std::string, Prof_var>;
        using Time_map = std::map<std::string, Time_var>;

        // Structure for columns.
        struct Column_struct
        {
            std::vector<int> coord;
            std::unique_ptr<Netcdf_file> data_file;
            std::unique_ptr<Netcdf_variable<int>> iter_var;
            std::unique_ptr<Netcdf_variable<TF>> time_var;
            Prof_map profs;
            Time_map time_series;
        };

        std::vector<Column_struct> columns;

        #ifdef USECUDA
        int* col_i_g;
        int* col_j_g;
        #endif

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
