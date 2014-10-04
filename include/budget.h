/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

#ifndef BUDGET
#define BUDGET

#include <string>

// forward declarations to reduce compilation time
class Model;
class Master;
class Input;
class Stats;
class Grid;
class Fields;
struct mask;

class Budget
{
  public:
    Budget(Model *, Input *);
    ~Budget();

    void init();
    void create();

    int execStats(mask *);

  private:
    Model  *model;
    Master *master;
    Stats  *stats;
    Grid   *grid;
    Fields *fields;

    std::string swbudget;

    double *umodel, *vmodel;

    int calcke(double *, double *, double *,
               double *, double *,
               double, double,
               double *, double *);

    int calctkebudget(double *, double *, double *, double *,
                      double *, double *,
                      double *, double *,
                      double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *,
                      double *, double *, double *,
                      double *, double *, double);

    int calctkebudget_buoy(double *, double *, double *, double *);

    int calcpe(double *, double *, double *, double *,
               double *,
               double *,
               double *, double *, double *,
               double *);

    int calcpebudget(double *, double *, double *, double *,
                     double *, double *, double *,
                     double *, double *, double *, double *,
                     double);

    int calcbpebudget(double *, double *, double *, double *, double *,
                      double *, double *, double *,
                      double *,
                      double *, double *, double *,
                      double);

    double calczsort   (double, double *, double *, int);
    double calcdzstardb(double, double *, double *);
};
#endif
