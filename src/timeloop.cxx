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

#include <iostream>
#include <cstdio>
#include <cmath>
#include "input.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "defines.h"
#include "constants.h"

template<typename TF>
Timeloop<TF>::Timeloop(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Input& input, const Sim_mode sim_mode) :
    master(masterin),
    grid(gridin),
    fields(fieldsin),
    ifactor(1e9)
{
    substep = 0;

    // obligatory parameters
    if (sim_mode == Sim_mode::Init)
    {
        starttime = 0.;
        phystarttime = input.get_item<double>("time", "phystarttime"  , "", 0.);
   }
    else
    {
        starttime = input.get_item<double>("time", "starttime", "");
        phystarttime = starttime + input.get_item<double>("time", "phystarttime"  , "", 0.);
    }

    datetime.tm_sec  = phystarttime;
    datetime.tm_year = 0; //default is 1900
    datetime.tm_mday = input.get_item<int>("time", "jday"  , "", 1);
    datetime.tm_isdst = -1;
    mktime ( &datetime );

    endtime  = input.get_item<double>("time", "endtime" , "");
    savetime = input.get_item<double>("time", "savetime", "");

    // optional parameters
    adaptivestep = input.get_item<bool>  ("time", "adaptivestep", "", true           );
    dtmax        = input.get_item<double>("time", "dtmax"       , "", Constants::dbig);
    dt           = input.get_item<double>("time", "dt"          , "", dtmax          );
    rkorder      = input.get_item<int>   ("time", "rkorder"     , "", 3              );
    outputiter   = input.get_item<int>   ("time", "outputiter"  , "", 20             );
    iotimeprec   = input.get_item<int>   ("time", "iotimeprec"  , "", 0              );

    if (sim_mode == Sim_mode::Post)
        postproctime = input.get_item<double>("time", "postproctime", "");

    // 3 and 4 are the only valid values for the rkorder
    if (!(rkorder == 3 || rkorder == 4))
    {
        master.print_error("\"%d\" is an illegal value for rkorder\n", rkorder);
        throw std::runtime_error("Illegal RK order value chosen");
    }

    // initializations
    loop      = true;
    time      = 0.;
    iteration = 0;

    // set or calculate all the integer times
    itime         = (unsigned long) 0;

    // add 0.5 to prevent roundoff errors
    iendtime      = (unsigned long)(ifactor * endtime + 0.5);
    istarttime    = (unsigned long)(ifactor * starttime + 0.5);
    idt           = (unsigned long)(ifactor * dt + 0.5);
    idtmax        = (unsigned long)(ifactor * dtmax + 0.5);
    isavetime     = (unsigned long)(ifactor * savetime + 0.5);
    if (sim_mode == Sim_mode::Post)
        ipostproctime = (unsigned long)(ifactor * postproctime + 0.5);

    idtlim = idt;

    // take the proper precision for the output files into account
    iiotimeprec = (unsigned long)(ifactor * std::pow(10., iotimeprec) + 0.5);

    // check whether starttime and savetime are an exact multiple of iotimeprec
    if ((istarttime % iiotimeprec) || (isavetime % iiotimeprec))
    {
        master.print_error("starttime or savetime is not an exact multiple of iotimeprec\n");
        throw 1;
    }

    iotime = (int)(istarttime / iiotimeprec);

    gettimeofday(&start, NULL);

    if (sim_mode == Sim_mode::Init)
        input.flag_as_used("time", "starttime", "");
}

template<typename TF>
Timeloop<TF>::~Timeloop()
{
}

template<typename TF>
void Timeloop<TF>::set_time_step_limit()
{
    idtlim = idtmax;

    // Check whether the run should be stopped because of the wall clock limit
    if (master.at_wall_clock_limit())
    {
        // Set the time step to the nearest multiple of iotimeprec
        idtlim = std::min(idtlim, iiotimeprec - itime % iiotimeprec);
    }

    idtlim = std::min(idtlim, isavetime - itime % isavetime);
}

template<typename TF>
void Timeloop<TF>::set_time_step_limit(unsigned long idtlimin)
{
    idtlim = std::min(idtlim, idtlimin);
}

template<typename TF>
void Timeloop<TF>::step_time()
{
    // Only step forward in time if we are not in a substep
    if (in_substep())
        return;

    time  += dt;
    itime += idt;
    iotime = (int)(itime/iiotimeprec);
    datetime.tm_sec = int(time + starttime);

    ++iteration;

    if (itime >= iendtime)
        loop = false;
}

template<typename TF>
bool Timeloop<TF>::do_check()
{
    if ((iteration % outputiter == 0 && !in_substep()) | !loop)
        return true;

    return false;
}

template<typename TF>
bool Timeloop<TF>::do_save()
{
    // Check whether the simulation has to stop due to the wallclock limit,
    // but only at a time step where actual saves can be made.
    if (itime % iiotimeprec == 0 && !in_substep() && master.at_wall_clock_limit())
    {
        master.print_warning("Simulation will be stopped after saving the restart files due to wall clock limit\n");

        // Stop looping
        loop = false;
        return true;
    }

    // Do not save directly after the start of the simulation and not in a substep
    if (itime % isavetime == 0 && iteration != 0 && !in_substep())
        return true;

    return false;
}

template<typename TF>
bool Timeloop<TF>::is_finished()
{
    // Return true if loop is false and vice versa.
    return !loop;
}

template<typename TF>
double Timeloop<TF>::check()
{
    gettimeofday(&end, NULL);

    double timeelapsed = (double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec-start.tv_usec) * 1.e-6;
    start = end;

    return timeelapsed;
}

template<typename TF>
void Timeloop<TF>::set_time_step()
{
    // Only set the time step if we are not in a substep
    if (in_substep())
        return;

    if (adaptivestep)
    {
        if (idt == 0)
        {
            master.print_error("Required time step less than precision %E of the time stepping\n", 1./ifactor);
            throw 1;
        }
        idt = idtlim;
        dt  = (double)idt / ifactor;
    }
}


namespace
{
    template<typename TF>
    void rk3(TF* restrict const a, TF* restrict const at, const int substep, const TF dt,
             const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
             const int jj, const int kk)
    {
        constexpr TF cA [] = {0., -5./9., -153./128.};
        constexpr TF cB [] = {1./3., 15./16., 8./15.};

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    a[ijk] += cB[substep]*dt*at[ijk];
                }

        const int substepn = (substep+1) % 3;

        // substep 0 resets the tendencies, because cA[0] == 0
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    at[ijk] *= cA[substepn];
                }
    }

    template<typename TF>
    void rk4(TF* restrict const a, TF* restrict const at, const int substep, const TF dt,
             const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
             const int jj, const int kk)
    {
        constexpr TF cA [] = {
            0.,
            - 567301805773./1357537059087.,
            -2404267990393./2016746695238.,
            -3550918686646./2091501179385.,
            -1275806237668./ 842570457699.};

        constexpr TF cB [] = {
            1432997174477./ 9575080441755.,
            5161836677717./13612068292357.,
            1720146321549./ 2090206949498.,
            3134564353537./ 4481467310338.,
            2277821191437./14882151754819.};

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];
                }

        const int substepn = (substep+1) % 5;

        // substep 0 resets the tendencies, because cA[0] == 0
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    at[ijk] = cA[substepn]*at[ijk];
                }
    }

    template<typename TF>
    inline TF rk3subdt(const TF dt, const int substep)
    {
        constexpr TF cB [] = {1./3., 15./16., 8./15.};
        return cB[substep]*dt;
    }

    template<typename TF>
    inline TF rk4subdt(const TF dt, const int substep)
    {
        constexpr TF cB [] = {
            1432997174477./ 9575080441755.,
            5161836677717./13612068292357.,
            1720146321549./ 2090206949498.,
            3134564353537./ 4481467310338.,
            2277821191437./14882151754819.};
        return cB[substep]*dt;
    }
}

#ifndef USECUDA
template<typename TF>
void Timeloop<TF>::exec()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (rkorder == 3)
    {
        for (auto& f : fields.at)
            rk3<TF>(fields.ap[f.first]->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

        substep = (substep+1) % 3;
    }

    if (rkorder == 4)
    {
        for (auto& f : fields.at)
            rk4<TF>(fields.ap[f.first]->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

        substep = (substep+1) % 5;
    }
}
#endif

template<typename TF>
double Timeloop<TF>::get_sub_time_step()
{
    // Value rkorder is 3 or 4, because it is checked in the constructor.
    if (rkorder == 3)
        return rk3subdt(dt, substep);
    else
        return rk4subdt(dt, substep);
}

template<typename TF>
bool Timeloop<TF>::in_substep()
{
    if (substep > 0)
        return true;
    else
        return false;
}

template<typename TF>
bool Timeloop<TF>::is_stats_step()
{
    // In case we are not in a substep and not at the first iteration
    // after a restart, we can could do statistics.
    if (!in_substep() && !((iteration > 0) && (itime == istarttime)))
        return true;
    else
        return false;
}

template<typename TF>
void Timeloop<TF>::save(int starttime)
{
    int nerror = 0;

    if (master.get_mpiid() == 0)
    {
        char filename[256];
        std::sprintf(filename, "time.%07d", starttime);

        master.print_message("Saving \"%s\" ... ", filename);

        FILE *pFile;
        pFile = fopen(filename, "wbx");

        if (pFile == NULL)
        {
            master.print_message("FAILED\n", filename);
            ++nerror;
        }
        else
        {
            fwrite(&itime    , sizeof(unsigned long), 1, pFile);
            fwrite(&idt      , sizeof(unsigned long), 1, pFile);
            fwrite(&iteration, sizeof(int), 1, pFile);

            fclose(pFile);
            master.print_message("OK\n");
        }
    }

    // Broadcast the error code to prevent deadlocks in case of error.
    master.broadcast(&nerror, 1);
    if (nerror)
        throw 1;
}

template<typename TF>
void Timeloop<TF>::load(int starttime)
{
    int nerror = 0;

    if (master.get_mpiid() == 0)
    {
        char filename[256];
        std::sprintf(filename, "time.%07d", starttime);

        master.print_message("Loading \"%s\" ... ", filename);

        FILE *pFile;
        pFile = fopen(filename, "rb");

        if (pFile == NULL)
        {
            master.print_error("\"%s\" does not exist\n", filename);
            ++nerror;
        }
        else
        {
            fread(&itime    , sizeof(unsigned long), 1, pFile);
            fread(&idt      , sizeof(unsigned long), 1, pFile);
            fread(&iteration, sizeof(int), 1, pFile);

            fclose(pFile);
        }
        master.print_message("OK\n");
    }

    master.broadcast(&nerror, 1);
    if (nerror)
        throw 1;

    master.broadcast(&itime    , 1);
    master.broadcast(&idt      , 1);
    master.broadcast(&iteration, 1);

    // calculate the double precision time from the integer time
    time = (double)itime / ifactor;
    dt   = (double)idt   / ifactor;
}

template<typename TF>
void Timeloop<TF>::step_post_proc_time()
{
    itime += ipostproctime;
    iotime = (int)(itime/iiotimeprec);

    if (itime > iendtime)
        loop = false;
}

template<typename TF>
Interpolation_factors<TF> Timeloop<TF>::get_interpolation_factors(const std::vector<double>& timevec)
{
    // 1. Get the indexes and factors for the interpolation in time
    std::vector<unsigned long> itimevec(timevec.size());
    for (int t=0; t<timevec.size(); ++t)
        itimevec[t] = static_cast<unsigned long>(ifactor * timevec[t] + 0.5);

    Interpolation_factors<TF> ifac;
    ifac.index1 = 0;
    for (auto& t : itimevec)
    {
        if (itime < t)
            break;
        else
            ++ifac.index1;
    }

    if (itime == itimevec[itimevec.size()-1])
        ifac.index1 = itimevec.size()-1;

    // 2. Calculate the weighting factor, accounting for out of range situations where the simulation is longer than the time range in input
    if (ifac.index1 == 0)
    {
        master.print_error("Interpolation time is out of range; t0 = %g; current time is %g \n", timevec[0], time);
        throw std::runtime_error("Interpolation time out of range");
    }
    else if (ifac.index1 == timevec.size())
    {
        master.print_error("Interpolation time is out of range; t1 = %g; current time is %g \n", timevec[timevec.size()-1], time);
        throw std::runtime_error("Interpolation time out of range");
    }
    else
    {
        ifac.index0 = ifac.index1-1;
        ifac.fac0 = TF(itimevec[ifac.index1] - itime) / TF(itimevec[ifac.index1] - itimevec[ifac.index0]);
        ifac.fac1 = TF(itime - itimevec[ifac.index0]) / TF(itimevec[ifac.index1] - itimevec[ifac.index0]);
    }

    return ifac;
}

template class Timeloop<double>;
template class Timeloop<float>;
