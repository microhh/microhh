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

#include <netcdf>
using namespace netCDF;

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;

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
        void set_mask_thres_pert(std::string, Field3d<TF>&, Field3d<TF>&, TF, Stats_mask_type );

        void exec(int, double, unsigned long);

        // Interface functions.
        void add_mask(const std::string);
        void add_prof(std::string, std::string, std::string, std::string);

        void add_fixed_prof(std::string, std::string, std::string, std::string, TF*);
        void add_time_series(std::string, std::string, std::string);

        void calc_stats(const std::string, const Field3d<TF>&, const int[3], const TF, const TF, const std::vector<std::string>);
        void set_prof(const std::string, const std::vector<TF>);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool swstats;           ///< Statistics on/off switch


        int statistics_counter;
        double sampletime;
        unsigned long isampletime;

        // Container for all stats, masks as uppermost in hierarchy
        Mask_map<TF> masks;
        std::vector<std::string> masklist;
        std::vector<unsigned int> mfield;
        std::vector<unsigned int> mfield_bot;

        void sanatize_operations_vector(std::vector<std::string>);
        bool wmean_set;

};
#endif
