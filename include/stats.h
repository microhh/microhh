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

#include <netcdfcpp.h>

class Master;
class Model;
class Grid;
class Fields;

// struct for profiles
struct ProfVar
{
  NcVar *ncvar;
  double *data;
};

// struct for time series
struct TimeSeriesVar
{
  NcVar *ncvar;
  double data;
};

// typedefs for containers of profiles and time series
typedef std::map<std::string, ProfVar> ProfMap;
typedef std::map<std::string, TimeSeriesVar> TimeSeriesMap;

// structure
struct Mask
{
  std::string name;
  NcFile *dataFile;
  NcDim  *z_dim, *zh_dim, *t_dim;
  NcVar  *t_var, *iter_var;
  ProfMap profs;
  TimeSeriesMap tseries;
};

typedef std::map<std::string, Mask> MaskMap;

class Stats
{
  public:
    Stats(Model *, Input *);
    ~Stats();

    void init(double);
    void create(int);

    unsigned long get_time_limit(unsigned long);
    void get_mask(Field3d *, Field3d *, Mask *);
    void exec(int, double, unsigned long);
    bool doStats();
    std::string getSwitch();

    // Container for all stats, masks as uppermost in hierarchy
    MaskMap masks;
    int *nmask;
    int *nmaskh;
    int nmaskbot;

    // Interface functions.
    void addMask(const std::string);
    void addProf(std::string, std::string, std::string, std::string);
    void addFixedProf(std::string, std::string, std::string, std::string, double *);
    void addTimeSeries(std::string, std::string, std::string);

    void calcArea   (double *, const int[3], int *);

    void calcMean  (double * const, const double * const,
                    const double, const int[3],
                    const double * const, const int * const);
    void calcMean2d(double * const, const double * const,
                    const double,
                    const double * const, const int * const);

    void calcMoment  (double *, double *, double *, double, const int[3], double *, int *);

    void calcDiff_2nd(double *, double *, double *, double, const int[3], double *, int *);
    void calcDiff_2nd(double *, double *, double *, double *, double *, double *, double *, double, const int[3], double *, int *);
    void calcDiff_4th(double *, double *, double *, double, const int[3], double *, int *);

    void calcGrad_2nd(double *, double *, double *, const int[3], double *, int *);
    void calcGrad_4th(double *, double *, double *, const int[3], double *, int *);

    void calcFlux_2nd(double *, double *, double *, double *, double *, double *, const int[3], double *, int *);
    void calcFlux_4th(double *, double *, double *, double *, const int[3], double *, int *);

    void addFluxes   (double *, double *, double *);
    void calcCount   (double *, double *, double, double *, int *);
    void calcPath    (double *, double *, int *, double *);
    void calcCover   (double *, double *, int *, double *, double);

    void calcSortedProf(double *, double *, double *);

  private:
    int nstats;

    // mask calculations
    void calcMask(double *, double *, double *, int *, int *, int *);

  protected:
    Model  *model;
    Grid   *grid;
    Fields *fields;
    Master *master;

    double sampletime;
    unsigned long isampletime;

    std::string swstats;

    static const int nthres = 0;
};
#endif
