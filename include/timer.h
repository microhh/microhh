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
        Timer(Master& masterin, const std::string& namein) : master(masterin), name(namein){};
        ~Timer(){};

        void add_timing(const std::string& name) {
            timings.emplace(name, Timing<TF>());
        };

        void start(const std::string& name) {
            timings.at(name).start_time = Time::now();
        };

        void stop(const std::string& name)
        {
            timings.at(name).end_time = Time::now();
            timings.at(name).elapsed +=
                (timings.at(name).end_time - timings.at(name).start_time).count()*TF(1e-9);
        };

        void save()
        {
            auto md = master.get_MPI_data();

            std::cout << "Timer \"" << name << "\"" << std::endl;
            for (auto& timer : timings)
            {
                // Get statistics along all MPI tasks
                TF min  = timer.second.elapsed;
                TF mean = timer.second.elapsed;
                TF max  = timer.second.elapsed;

                master.sum(&mean, 1);
                master.min(&min,  1);
                master.max(&max,  1);

                mean /= md.nprocs;

                if (md.mpiid == 0)
                    std::cout << " - " << timer.first << " min=" << min << " mean=" << mean << " max=" << max << std::endl;

                // Reset timer.
                timer.second.elapsed = 0;
            }
        };

    private:
        Master& master;

        std::string name;
        std::map<std::string, Timing<TF>> timings;
};
#endif
