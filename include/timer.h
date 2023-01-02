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

#ifndef TIMER_H
#define TIMER_H

#include <map>
#include <string>
#include <chrono>
#include <iostream>

#include "netcdf_interface.h"

#ifdef USEMPI
#include <mpi.h>
#endif

class Master;

namespace
{
    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::time_point<Time> time_point;
}

template<typename TF>
struct Timing
{
    time_point start_time;
    time_point end_time;
    TF elapsed = 0;
};


template<typename TF>
class Timer
{
    public:
        Timer(Master& masterin, Input& inputin, const std::string& namein) : master(masterin), name(namein)
        {
            sw_timer = inputin.get_item<bool>("timer", "swtimer", "", false);
        };

        ~Timer(){};

        void create()
        {
            if (sw_timer)
            {
                // Create NetCDF file.
                std::string file_name = this->name + "_timings.nc";

                nc_file = std::make_unique<Netcdf_file>(master, file_name, Netcdf_mode::Create);
                nc_file->add_dimension("time");

                time_var = std::make_unique<Netcdf_variable<TF>>(
                        nc_file->template add_variable<TF>("time", {"time"}));
                time_var->add_attribute("units", "seconds since start");
                time_var->add_attribute("long_name", "Time");

                nc_file->sync();
            }
        };

        void add_timing(const std::string& name)
        {
            if (sw_timer)
            {
                timings.emplace(name, Timing<TF>());

                // Add NetCDF variable.
                std::vector<std::string> modes = {"min", "mean", "max"};
                for (auto& mode : modes)
                {
                    std::string var_name = name + "_" + mode;

                    Time_var var{nc_file->template add_variable<TF>(var_name, {"time"}), TF(0)};
                    var.ncvar.add_attribute("units", "s");

                    time_series.emplace(
                            std::piecewise_construct, std::forward_as_tuple(var_name), std::forward_as_tuple(std::move(var)));
                }
                nc_file->sync();
            }
        };

        void start(const std::string& name)
        {
            if (sw_timer)
                timings.at(name).start_time = Time::now();
        };

        void stop(const std::string& name)
        {
            if (sw_timer)
            {
                timings.at(name).end_time=Time::now();
                timings.at(name).elapsed+=
                        (timings.at(name).end_time - timings.at(name).start_time).count() * TF(1e-9);
            }
        };

        void save(const int time)
        {
            if (sw_timer)
            {
                auto md = master.get_MPI_data();

                for (auto& timer : timings)
                {
                    // Get statistics along all MPI tasks
                    TF min = timer.second.elapsed;
                    TF mean = timer.second.elapsed;
                    TF max = timer.second.elapsed;

                    master.sum(&mean, 1);
                    master.min(&min, 1);
                    master.max(&max, 1);

                    mean/=md.nprocs;

                    std::string name_min  = timer.first + "_min";
                    std::string name_mean = timer.first + "_mean";
                    std::string name_max  = timer.first + "_max";

                    time_series.at(name_min) .data = min;
                    time_series.at(name_mean).data = mean;
                    time_series.at(name_max) .data = max;
                }

                // Store value in NetCDF file.
                const std::vector<int> time_index{statistics_counter};
                time_var->insert(time, time_index);

                for (auto& tser : time_series)
                    tser.second.ncvar.insert(tser.second.data, time_index);

                nc_file->sync();
                ++statistics_counter;

                // Reset timers.
                for (auto& timer : timings)
                    timer.second.elapsed = 0;
            }
        };

    private:
        Master& master;

        bool sw_timer;
        std::string name;
        std::map<std::string, Timing<TF>> timings;

        // NetCDF output
        std::unique_ptr<Netcdf_file> nc_file;
        std::unique_ptr<Netcdf_variable<TF>> time_var;
        int statistics_counter = 0;

        struct Time_var
        {
            Netcdf_variable<TF> ncvar;
            TF data;
        };
        std::map<std::string, Time_var> time_series;
};
#endif
