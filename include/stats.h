/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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

//#include <netcdfcpp.h>
#include <netcdf>
using namespace netCDF;

class Master;
class Model;
class Grid;
class Fields;

// struct for profiles
struct Prof_var
{
    NcVar ncvar;
    double* data;
};

// struct for time series
struct Time_series_var
{
    NcVar ncvar;
    double data;
};

// typedefs for containers of profiles and time series
typedef std::map<std::string, Prof_var> Prof_map;
typedef std::map<std::string, Time_series_var> Time_series_map;

// structure
struct Mask
{
    std::string name;
    NcFile* dataFile;
    NcDim z_dim;
    NcDim zh_dim;
    NcDim t_dim;
    NcVar iter_var;
    NcVar t_var;
    Prof_map profs;
    Time_series_map tseries;
};

typedef std::map<std::string, Mask> Mask_map;

class Stats
{
    public:
        Stats(Model*, Input*);
        ~Stats();

        void init(double);
        void create(int);

        unsigned long get_time_limit(unsigned long);
        void get_mask(Field3d*, Field3d*, Mask*);
        void exec(int, double, unsigned long);
        bool doStats();
        std::string get_switch();

        // Container for all stats, masks as uppermost in hierarchy
        Mask_map masks;
        int* nmask;
        int* nmaskh;
        int nmaskbot;

        // Interface functions.
        void add_mask(const std::string);
        void add_prof(std::string, std::string, std::string, std::string);
        void add_fixed_prof(std::string, std::string, std::string, std::string, double*);
        void add_time_series(std::string, std::string, std::string);

        void calc_area(double*, const int[3], int*);

        void calc_mean(double* const, const double* const,
                       const double, const int[3],
                       const double* const, const int* const);

        void calc_mean2d(double* const, const double* const,
                         const double,
                         const double* const, const int* const);

        void calc_moment  (double*, double*, double*, double, const int[3], double*, int*);

        void calc_diff_2nd(double*, double*, double*, double, const int[3], double*, int*);
        void calc_diff_2nd(double*, double*, double*, double*, double*,
                           double*, double*, double, const int[3], double*, int*);
        void calc_diff_4th(double*, double*, double*, double, const int[3], double*, int*);

        void calc_grad_2nd(double*, double*, double*, const int[3], double*, int*);
        void calc_grad_4th(double*, double*, double*, const int[3], double*, int*);

        void calc_flux_2nd(double*, double*, double*, double*, double*, double*, const int[3], double*, int*);
        void calc_flux_4th(double*, double*, double*, double*, const int[3], double*, int*);

        void add_fluxes   (double*, double*, double*);
        void calc_count   (double*, double*, double, double*, int*);
        void calc_path    (double*, double*, int*, double*);
        void calc_cover   (double*, double*, int*, double*, double);

        void calc_sorted_prof(double*, double*, double*);

    private:
        int nstats;

        // mask calculations
        void calc_mask(double*, double*, double*, int*, int*, int*);

    protected:
        Model*  model;
        Grid*   grid;
        Fields* fields;
        Master* master;

        double sampletime;
        unsigned long isampletime;

        std::string swstats;

        static const int nthres = 0;
};
#endif
