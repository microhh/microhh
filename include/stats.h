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

#ifndef STATS
#define STATS

#include <regex>
#include <netcdf>
using namespace netCDF;

#include "boundary_cyclic.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Diff;

// Struct for profiles
template<typename TF>
struct Prof_var
{
    NcVar ncvar;
    std::vector<TF> data;
};

// Struct for time series
template<typename TF>
struct Time_series_var
{
    NcVar ncvar;
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
    NcFile* data_file;
    NcDim z_dim;
    NcDim zh_dim;
    NcDim t_dim;
    NcVar iter_var;
    NcVar t_var;
    Prof_map<TF> profs;
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
        Stats(Master&, Grid<TF>&, Fields<TF>&, Input&);  ///< Constructor of the statistics class
        ~Stats();

        void init(double);
        void create(int, std::string);

        unsigned long get_time_limit(unsigned long);
        bool get_switch() { return swstats; }
        bool do_statistics(unsigned long);

        void initialize_masks();
        void finalize_masks();

        const std::vector<std::string>& get_mask_list();
        void set_mask_thres(std::string, Field3d<TF>&, Field3d<TF>&, TF, Stats_mask_type );

        void exec(int, double, unsigned long);

        // Interface functions.
        void add_mask(const std::string);
        void add_prof(std::string, std::string, std::string, std::string, Stats_whitelist_type = Stats_whitelist_type::Default);

        void add_fixed_prof(std::string, std::string, std::string, std::string, TF*);
        void add_time_series(std::string, std::string, std::string, Stats_whitelist_type = Stats_whitelist_type::Default);

        void calc_stats(const std::string, const Field3d<TF>&, const TF, const TF, std::vector<std::string>, Diff<TF>&);
        void calc_stats_mean(const std::string, const Field3d<TF>&, const TF, const TF);
        void calc_stats_2d(const std::string, const std::vector<TF>&, const TF, std::vector<std::string>);
        void calc_covariance(const std::string, const Field3d<TF>&, const TF, const TF, const int,
                             const std::string, const Field3d<TF>&, const TF, const TF, const int);
        void set_prof(const std::string, const std::vector<TF>);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Boundary_cyclic<TF> boundary_cyclic;

        bool swstats;           ///< Statistics on/off switch
        std::vector<std::regex> whitelist;
        std::vector<std::regex> blacklist;
        std::vector<std::string> varlist;
        bool is_blacklisted(std::string, Stats_whitelist_type);

        int statistics_counter;
        double sampletime;
        unsigned long isampletime;

        // Container for all stats, masks as uppermost in hierarchy
        Mask_map<TF> masks;
        std::vector<std::string> masklist;
        std::vector<unsigned int> mfield;
        std::vector<unsigned int> mfield_bot;

        void calc_flux_2nd(TF*, const TF* const, const TF* const, TF* const, const TF* const, TF*, const int*, const unsigned int* const, const unsigned int, const int* const,
                          const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_flux_4th(TF*, const TF* const, const TF* const, TF* const, const TF* const, TF*, const int*, const unsigned int* const, const unsigned int, const int* const,
                        const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_grad_2nd(TF* const restrict, const TF* const restrict, const TF* const restrict,
                            const unsigned int* const, const unsigned int, const int* const,
                            const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_grad_4th(TF* const restrict, const TF* const restrict, const TF* const restrict,
                const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_diff_2nd(TF* restrict, TF* restrict, const TF* restrict, TF, const int*,
            const unsigned int* const, const unsigned int, const int* const,
            const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_diff_2nd(TF* restrict, TF* restrict w, TF* restrict,
                TF* restrict, const TF* restrict,
                TF* restrict, TF* restrict, const TF, const int*,
                const unsigned int* const, const unsigned int, const int* const,
                const TF, const TF, const int, const int, const int, const int, const int, const int, const int, const int);
        void calc_diff_4th(
                TF* restrict, TF* restrict, const TF* restrict,
                const TF, const int*,
                const unsigned int* const, const unsigned int, const int* const,
                const int, const int, const int, const int, const int, const int, const int, const int);

        void sanitize_operations_vector(std::vector<std::string>);
        bool wmean_set;

};
#endif
