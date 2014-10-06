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

#ifndef TIMELOOP
#define TIMELOOP

#include <sys/time.h>
#include <string>

class Input;
class Master;
class Model;
class Grid;
class Fields;

class Timeloop
{
  public:
    Timeloop(Model *, Input *);
    ~Timeloop();

    void stepTime();
    int postprocstep();
    void setTimeStep();

    int exec();

    // Query functions for main loop
    bool inSubStep();
    bool isStatsStep();
    bool doCheck();
    bool doSave();

    double getSubTimeStep();

    int settimelim();

    unsigned long settimelim(unsigned long);

    void save(int);
    void load(int);

    double check();


    // variables
    int substep;
    bool loop;
    bool adaptivestep;

    double time;
    double dt;
    double dtmax;
    double endtime;
    double postproctime;
    double savetime;
    double starttime;
    double dtlim;

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

    double ifactor;

  private:
    Master *master;
    Model  *model;
    Grid   *grid;
    Fields *fields;

    timeval start;
    timeval end;

    int rkorder;

    int outputiter;

    int rk3(double *, double *, double);
    int rk4(double *, double *, double);
    double rk3subdt(double);
    double rk4subdt(double);

    int rk3_GPU(double *, double *, double);
    int rk4_GPU(double *, double *, double);
};
#endif
