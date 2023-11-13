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

#ifndef STATS_H
#define STATS_H

#include <regex>
#include "boundary_cyclic.h"

class Master;
class Input;
class Netcdf_file;
template<typename> class Grid;
template<typename> class Soil_grid;
template<typename> class Background;
template<typename> class Fields;
template<typename> class Field3d;
template<typename> class Advec;
template<typename> class Diff;
template<typename> class Timeloop;
template<typename> class Netcdf_variable;
template<typename> class Field3d;

// Struct for profiles
enum class Level_type {Full, Half};

template<typename TF>
struct Prof_var
{
    Netcdf_variable<TF> ncvar;
    std::vector<TF> data;
    Level_type level;
};

// Struct for time series
template<typename TF>
struct Time_series_var
{
    Netcdf_variable<TF> ncvar;
    TF data;
};

// Typedefs for containers of profiles and time series
template<typename TF>
using Prof_map = std::map<std::string, Prof_var<TF>>;

template<typename TF>
using Time_series_map = std::map<std::string, Time_series_var<TF>>;

// structure
template<typename TF>
struct Mask
{
    std::string name;
    unsigned int flag;
    unsigned int flagh;
    std::vector<int> nmask;
    std::vector<int> nmaskh;
    int nmask_bot;

    std::unique_ptr<Netcdf_file> data_file;
    std::unique_ptr<Netcdf_variable<int>> iter_var;
    std::unique_ptr<Netcdf_variable<TF>> time_var;
    Prof_map<TF> profs;
    Prof_map<TF> soil_profs;
    Prof_map<TF> background_profs;
    Time_series_map<TF> tseries;
};

template<typename TF>
using Mask_map = std::map<std::string, Mask<TF>>;

enum class Stats_mask_type {Plus, Min};
enum class Stats_whitelist_type {White, Black, Default};

template<typename TF>
class Stats
{
    public:
        Stats(Master&, Grid<TF>&, Soil_grid<TF>&, Background<TF>&, Fields<TF>&, Advec<TF>&, Diff<TF>&, Input&);
        ~Stats();

        void init(double);
        void create(const Timeloop<TF>&, std::string);

        unsigned long get_time_limit(unsigned long);
        bool get_switch() { return swstats; }
        bool do_statistics(unsigned long);
        bool do_tendency() {return swtendency; }
        void set_tendency(bool);

        void initialize_masks();
        void finalize_masks();

        const std::vector<std::string>& get_mask_list();
        void set_mask_thres(std::string, Field3d<TF>&, Field3d<TF>&, TF, Stats_mask_type );

        void exec(const int, const double, const unsigned long);

        // Interface functions.
        void add_dimension(const std::string&, const int);
        void add_mask(const std::string&);
        void add_prof(
                const std::string&, const std::string&,
                const std::string&, const std::string&, const std::string&,
                Stats_whitelist_type = Stats_whitelist_type::Default);
        void add_profs(
                const Field3d<TF>&, const std::string&, std::vector<std::string>, const std::string&);
        void add_tendency(
                const Field3d<TF>&, const std::string&, const std::string&,
                const std::string&, const std::string& group_name="tend");

        void add_covariance(
                const Field3d<TF>&, const Field3d<TF>&, const std::string&, const std::string& group_name="covariances");

        void add_fixed_prof(
                const std::string&, const std::string&,
                const std::string&, const std::string&, const std::string&,
                const std::vector<TF>&);

        void add_fixed_prof_raw(
                const std::string&, const std::string&,
                const std::string&, const std::string&, const std::string&,
                const std::vector<TF>&);

        void add_time_series(
                const std::string&, const std::string&,
                const std::string&, const std::string&, Stats_whitelist_type=Stats_whitelist_type::Default);

        void calc_mask_stats(
                std::pair<const std::string, Mask<TF>>&,
                const std::string&, const Field3d<TF>&, const TF, const TF);

        void calc_mask_mean_profile(
                std::vector<TF>&,
                const std::pair<const std::string, Mask<TF>>&,
                const Field3d<TF>&);

        void calc_stats(const std::string&, const Field3d<TF>&, const TF, const TF);
        void calc_stats_mean(const std::string&, const Field3d<TF>&, const TF);
        void calc_stats_moments(const std::string&, const Field3d<TF>&, const TF);
        void calc_stats_w(const std::string&, const Field3d<TF>&, const TF);
        void calc_stats_diff(const std::string&, const Field3d<TF>&, const TF);
        void calc_stats_flux(const std::string&, const Field3d<TF>&, const TF);
        void calc_stats_grad(const std::string&, const Field3d<TF>&);
        void calc_stats_path(const std::string&, const Field3d<TF>&);
        void calc_stats_cover(const std::string&, const Field3d<TF>&, const TF, const TF);
        void calc_stats_frac(const std::string&, const Field3d<TF>&, const TF, const TF);

        void calc_stats_2d(const std::string&, const std::vector<TF>&, const TF);
        void calc_stats_soil(const std::string, const std::vector<TF>&, const TF);
        void calc_covariance(const std::string&, const Field3d<TF>&, const TF, const TF, const int,
                             const std::string&, const Field3d<TF>&, const TF, const TF, const int);
        void calc_tend(Field3d<TF>&, const std::string&);
        void set_prof(const std::string&, const std::vector<TF>&);
        void set_prof_background(const std::string&, const std::vector<TF>&);
        void set_time_series(const std::string&, const TF);

        Mask_map<TF>& get_masks() { return masks; }

    private:
        Master& master;
        Grid<TF>& grid;
        Soil_grid<TF>& soil_grid;
        Background<TF>& background;
        Fields<TF>& fields;
        Advec<TF>& advec;
        Diff<TF>& diff;
        Boundary_cyclic<TF> boundary_cyclic;

        bool swstats;           ///< Statistics on/off switch
        bool swtendency;
        bool doing_tendency;
        std::vector<std::regex> whitelist;
        std::vector<std::regex> blacklist;
        std::vector<std::string> varlist;
        std::vector<std::string> varlist_soil;
        std::vector<std::string> varlist_background;
        
        void add_operation(std::vector<std::string>&, const std::string&, const std::string&);
        void sanitize_operations_vector(const std::string&, std::vector<std::string>&);
        bool is_blacklisted(const std::string&, Stats_whitelist_type = Stats_whitelist_type::Default);

        int statistics_counter;
        double sampletime;
        unsigned long isampletime;

        // Container for all stats, masks as uppermost in hierarchy
        Mask_map<TF> masks;
        std::vector<std::string> masklist;
        std::vector<unsigned int> mfield;
        std::vector<unsigned int> mfield_bot;

        // Tendency calculations
        std::map<std::string, std::vector<std::string>> tendency_order;

        void calc_flux_2nd(
                TF*, const TF* const, const TF* const, const TF, TF* const, const TF* const, TF*,
                const int*, const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_flux_4th(
                TF*, const TF* const, const TF* const, TF* const, const TF* const, TF*,
                const int*, const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_grad_2nd(
                TF* const restrict, const TF* const restrict, const TF* const restrict,
                const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_grad_4th(
                TF* const restrict, const TF* const restrict, const TF* const restrict,
                const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_diff_2nd(
                TF* restrict, TF* restrict, const TF* restrict, TF, const int*,
                const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_diff_2nd(
                TF* restrict, TF* restrict w, TF* restrict,
                TF* restrict, const TF* restrict,
                TF* restrict, TF* restrict, const TF, const int*,
                const unsigned int* const, const unsigned int, const int* const,
                const TF, const TF, const int, const int, const int, const int,
                const int, const int, const int, const int);
        void calc_diff_4th(
                TF* restrict, TF* restrict, const TF* restrict,
                const TF, const int*,
                const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);

        bool wmean_set;

};
#endif
