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

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

class Master;
class Input;
class Netcdf_file;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;
template<typename> class Netcdf_variable;

template<typename TF>
struct Time_var
{
    Netcdf_variable<TF> ncvar;
    TF data;
};

template<typename TF>
struct Single_trajectory
{
    std::string name;    // Descriptive name of tile

    // Input vectors with time and (x,y,z) location.
    std::vector<double> time_in;
    std::vector<TF> x_in;
    std::vector<TF> y_in;
    std::vector<TF> z_in;

    // Current location of trajectory.
    TF x_loc;
    TF y_loc;
    TF z_loc;

    // NetCDF file.
    std::unique_ptr<Netcdf_file> data_file;

    // Individual time series of sampled statistics.
    std::map<std::string, Time_var<TF>> time_series;

    // NetCDF variables of time/location.
    std::unique_ptr<Netcdf_variable<TF>> time_var;
    std::unique_ptr<Netcdf_variable<TF>> x_var;
    std::unique_ptr<Netcdf_variable<TF>> y_var;
    std::unique_ptr<Netcdf_variable<TF>> z_var;
};

template<typename TF>
using Trajectory_map = std::map<std::string, Single_trajectory<TF>>;


template<typename TF>
class Trajectory
{
    public:
        Trajectory(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Trajectory();

        void init();
        void create(Input&, Netcdf_file&, Timeloop<TF>&, std::string);

        unsigned long get_time_limit(unsigned long);
        bool do_trajectory(unsigned long);
        bool get_switch() { return sw_trajectory; }
        void exec(Timeloop<TF>&, double, unsigned long);

    private:
        std::vector<std::string> names;     // Individual trajectories.
        std::vector<std::string> variables; // Names of prognostic/diagnostic variables to sample.

        Trajectory_map<TF> trajectories;


    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_trajectory;

        int statistics_counter;
        double sampletime;
        unsigned long isampletime;
};
#endif
