/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef TIMELOOP_H
#define TIMELOOP_H

#include <string>
#include <vector>
#include <ctime>

class Master;
template<typename> class Grid;
template<typename> class Soil_grid;
template<typename> class Fields;
class Input;
enum class Sim_mode;


template<typename TF>
struct Interpolation_factors
{
    unsigned int index0;
    unsigned int index1;
    TF fac0;
    TF fac1;
};

template<typename TF>
class Timeloop
{
    public:
        Timeloop(Master&, Grid<TF>&, Soil_grid<TF>&,
                Fields<TF>&, Input&, const Sim_mode);
        ~Timeloop();

        void step_time();
        void step_post_proc_time();
        void set_time_step();
        void set_time_step_limit();
        void set_time_step_limit(unsigned long);
        double get_sub_time_step() const;

        Interpolation_factors<TF> get_interpolation_factors(const std::vector<double>&);

        void exec();

        double check();

        void save(int, unsigned long, unsigned long, int);
        void load(int);

        // Query functions for main loop
        bool in_substep();
        bool is_stats_step();
        bool do_check();
        bool do_save();
        bool is_finished();

        // Accessors for other classes
        double get_time() const { return time; }
        double get_dt() const { return dt; }
        double get_ifactor() const { return ifactor; }
        unsigned long get_iiotimeprec() const { return iiotimeprec; }
        unsigned long get_itime() const { return itime; }
        unsigned long get_idt() const { return idt; }
        unsigned long get_isavetime() const { return isavetime; }
        int get_iotime() const { return iotime; }
        int get_iteration() const { return iteration; }
        int get_substep() const { return substep; }

        // Functions for UTC time support.
        bool has_utc_time() const { return flag_utc_time; }
        std::string get_datetime_utc_start_string() const;
        double calc_day_of_year() const;
        double calc_hour_of_day() const;
        int get_year() const;

    private:
        Master& master;
        Grid<TF>& grid;
        Soil_grid<TF>& soil_grid;
        Fields<TF>& fields;

        timeval start;
        timeval end;

        int rkorder;
        int outputiter;

        // Variables
        bool loop;

        int substep;
        bool adaptivestep;

        double dt;
        double dtmax;

        double time;
        double starttime;
        double endtime;
        double savetime;
        double postproctime;
        bool flag_utc_time;
        std::tm tm_utc_start;

        int iteration;
        int iotime;
        int iotimeprec;

        unsigned long itime;
        unsigned long istarttime;
        unsigned long iendtime;
        unsigned long idt;
        unsigned long idtmax;
        unsigned long ipostproctime;
        unsigned long isavetime;
        unsigned long idtlim;
        unsigned long iiotimeprec;

        const double ifactor;
};
#endif
