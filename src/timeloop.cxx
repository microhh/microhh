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

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <ctime>
#include <sys/time.h>

#include "input.h"
#include "master.h"
#include "grid.h"
#include "soil_grid.h"
#include "soil_field3d.h"
#include "fields.h"
#include "timeloop.h"
#include "defines.h"
#include "constants.h"

template<typename TF>
Timeloop<TF>::Timeloop(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Input& input, const Sim_mode sim_mode) :
    master(masterin),
    grid(gridin),
    soil_grid(soilgridin),
    fields(fieldsin),
    flag_utc_time(false),
    ifactor(1e9)
{
    setenv("TZ", "utc", 1);

    substep = 0;

    // Obligatory parameters.
    if (sim_mode == Sim_mode::Init)
        starttime = 0.;
    else
        starttime = input.get_item<double>("time", "starttime", "");

    endtime  = input.get_item<double>("time", "endtime" , "");
    savetime = input.get_item<double>("time", "savetime", "");

    // Optional parameters.
    adaptivestep = input.get_item<bool>  ("time", "adaptivestep", "", true           );
    dtmax        = input.get_item<double>("time", "dtmax"       , "", Constants::dbig);
    dt           = input.get_item<double>("time", "dt"          , "", dtmax          );
    rkorder      = input.get_item<int>   ("time", "rkorder"     , "", 3              );
    outputiter   = input.get_item<int>   ("time", "outputiter"  , "", 20             );
    iotimeprec   = input.get_item<int>   ("time", "iotimeprec"  , "", 0              );

    // Get a datetime in UTC.
    std::string datetime_utc_string = input.get_item<std::string>("time", "datetime_utc", "", "");
    if (datetime_utc_string != "")
    {
        flag_utc_time = true;
        strptime(datetime_utc_string.c_str(), "%Y-%m-%d %H:%M:%S", &tm_utc_start);

        // NOTE: the following fields are NOT set by `strptime()`, which can lead to undefined behaviour.
        tm_utc_start.tm_isdst = 0;      // no daylight saving offset.
        tm_utc_start.tm_gmtoff = 0;     // no offset from UTC.
        tm_utc_start.tm_zone = "utc";   // time zone = UTC.
    }

    if (sim_mode == Sim_mode::Post)
        postproctime = input.get_item<double>("time", "postproctime", "");

    // 3 and 4 are the only valid values for the rkorder
    if (!(rkorder == 3 || rkorder == 4))
    {
        std::string msg = std::to_string(rkorder) + " is an illegal value for rkorder";
        throw std::runtime_error(msg);
    }

    // initializations
    loop      = true;
    time      = 0.;
    iteration = 0;

    // set or calculate all the integer times
    itime      = static_cast<unsigned long>(0);

    // add 0.5 to prevent roundoff errors
    iendtime   = static_cast<unsigned long>(ifactor * endtime + 0.5);
    istarttime = static_cast<unsigned long>(ifactor * starttime + 0.5);
    idt        = static_cast<unsigned long>(ifactor * dt + 0.5);
    idtmax     = static_cast<unsigned long>(ifactor * dtmax + 0.5);
    isavetime  = static_cast<unsigned long>(ifactor * savetime + 0.5);
    if (sim_mode == Sim_mode::Post)
        ipostproctime = static_cast<unsigned long>(ifactor * postproctime + 0.5);

    idtlim = idt;

    // take the proper precision for the output files into account
    iiotimeprec = static_cast<unsigned long>(ifactor * std::pow(10., iotimeprec) + 0.5);

    // check whether starttime and savetime are an exact multiple of iotimeprec
    if ((istarttime % iiotimeprec) || (isavetime % iiotimeprec))
    {
        std::string msg = " Starttime or savetime is not an exact multiple of iotimeprec";
        throw std::runtime_error(msg);
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
    if (itime < iendtime)
        idtlim = std::min(idtlim, iendtime - itime);
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
    iotime = static_cast<int>(itime/iiotimeprec);

    ++iteration;

    if (itime >= iendtime)
        loop = false;
}

template<typename TF>
bool Timeloop<TF>::do_check()
{
    // Print every RK3 substep if outputiter == 0. Useful for debugging.
    if (outputiter == 0)
        return true;
    else if ((iteration % outputiter == 0 && !in_substep()) | !loop)
        return true;
    else
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

    double timeelapsed = static_cast<double>(end.tv_sec-start.tv_sec) + static_cast<double>(end.tv_usec-start.tv_usec) * 1.e-6;
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
            std::string msg = "Required time step less than precision " + std::to_string(1./ifactor) + " of the time stepping";
            throw std::runtime_error(msg);
        }
        idt = idtlim;
        dt  = static_cast<double>(idt) / ifactor;
    }
}

namespace
{
    template<typename TF>
    void rk3(TF* restrict const a, TF* restrict const at, const int substep, const TF dt,
             const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
             const int jj, const int kk, const int ncells)
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
        if (substepn == 0)
        {
            for (int n=0; n<ncells; ++n)
                at[n] = TF(0.);
        }
        else
        {
            for (int k=kstart; k<kend; ++k)
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        at[ijk] *= cA[substepn];
                    }
        }
    }

    template<typename TF>
    void rk4(TF* restrict const a, TF* restrict const at, const int substep, const TF dt,
             const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
             const int jj, const int kk, const int ncells)
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
        if (substepn == 0)
        {
            for (int n=0; n<ncells; ++n)
                at[n] = TF(0.);
        }
        else
        {
            for (int k=kstart; k<kend; ++k)
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        at[ijk] = cA[substepn]*at[ijk];
                    }
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
    auto& gd  = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int kstart_2d = 0;
    const int kend_2d = 1;

    if (rkorder == 3)
    {
        // Atmospheric fields
        for (auto& f : fields.at)
            rk3<TF>(fields.ap.at(f.first)->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells, gd.ncells);

        // Soil fields
        for (auto& f : fields.sts)
            rk3<TF>(fields.sps.at(f.first)->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, sgd.kstart, sgd.kend,
                    gd.icells, gd.ijcells, sgd.ncells);

        // 2D fields
        for (auto& f : fields.at2d)
            rk3<TF>(fields.ap2d.at(f.first)->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, kstart_2d, kend_2d,
                    gd.icells, gd.ijcells, gd.ijcells);

        substep = (substep+1) % 3;
    }

    if (rkorder == 4)
    {
        // Atmospheric fields
        for (auto& f : fields.at)
            rk4<TF>(fields.ap.at(f.first)->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells, gd.ncells);

        // Soil fields
        for (auto& f : fields.sts)
            rk4<TF>(fields.sps.at(f.first)->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, sgd.kstart, sgd.kend,
                    gd.icells, gd.ijcells, sgd.ncells);

        // 2D fields
        for (auto& f : fields.at2d)
            rk4<TF>(fields.ap2d.at(f.first)->fld.data(), f.second->fld.data(), substep, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, kstart_2d, kend_2d,
                    gd.icells, gd.ijcells, gd.ijcells);

        substep = (substep+1) % 5;
    }
}
#endif

template<typename TF>
double Timeloop<TF>::get_sub_time_step() const
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
void Timeloop<TF>::save(int starttime, unsigned long itime_in, unsigned long idt_in, int iteration_in)
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
            fwrite(&itime_in    , sizeof(unsigned long), 1, pFile);
            fwrite(&idt_in      , sizeof(unsigned long), 1, pFile);
            fwrite(&iteration_in, sizeof(int), 1, pFile);

            fclose(pFile);
            master.print_message("OK\n");
        }
    }

    // Broadcast the error code to prevent deadlocks in case of error.
    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("In Timeloop::save");
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
            std::string msg = std::string(filename) + " does not exist";
            master.print_message(msg);
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
    // Broadcast the error code to prevent deadlocks in case of error.
    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("In Timeloop::load");

    master.broadcast(&itime    , 1);
    master.broadcast(&idt      , 1);
    master.broadcast(&iteration, 1);

    // calculate the double precision time from the integer time
    time = static_cast<double>(itime) / ifactor;
    dt   = static_cast<double>(idt)   / ifactor;
}

template<typename TF>
void Timeloop<TF>::step_post_proc_time()
{
    itime += ipostproctime;
    iotime = static_cast<int>(itime/iiotimeprec);

    if (itime > iendtime)
        loop = false;
}

namespace
{
    std::tm calc_tm_actual(const std::tm& tm_start, const double time)
    {
        std::tm tm_actual = tm_start;
        tm_actual.tm_sec += static_cast<int>(time);
        std::mktime(&tm_actual);
        return tm_actual;
    }
}

template<typename TF>
double Timeloop<TF>::calc_hour_of_day() const
{
    if (!flag_utc_time)
        throw std::runtime_error("No datetime in UTC specified");

    std::tm tm_actual = calc_tm_actual(tm_utc_start, time);
    const double frac_hour = ( tm_actual.tm_min*60
                             + tm_actual.tm_sec + std::fmod(time, 1.) ) / 3600.;
    return tm_actual.tm_hour + frac_hour; // Counting starts at 0 in std::tm, thus add 1.
}

template<typename TF>
double Timeloop<TF>::calc_day_of_year() const
{
    if (!flag_utc_time)
        throw std::runtime_error("No datetime in UTC specified");

    std::tm tm_actual = calc_tm_actual(tm_utc_start, time);
    const double frac_day = ( tm_actual.tm_hour*3600.
                            + tm_actual.tm_min*60.
                            + tm_actual.tm_sec + std::fmod(time, 1.) ) / 86400.;
    return tm_actual.tm_yday+1. + frac_day; // Counting starts at 0 in std::tm, thus add 1.
}

template<typename TF>
int Timeloop<TF>::get_year() const
{
    if (!flag_utc_time)
        throw std::runtime_error("No datetime in UTC specified");

    std::tm tm_actual = calc_tm_actual(tm_utc_start, time);
    return tm_actual.tm_year;
}

template<typename TF>
std::string Timeloop<TF>::get_datetime_utc_start_string() const
{
    if (!flag_utc_time)
        throw std::runtime_error("No datetime in UTC specified");

    std::ostringstream ss;

    // Year is relative to 1900, month count starts at 0.
    ss << std::setfill('0') << std::setw(4) << tm_utc_start.tm_year+1900 << "-";
    ss << std::setfill('0') << std::setw(2) << tm_utc_start.tm_mon+1     << "-";
    ss << std::setfill('0') << std::setw(2) << tm_utc_start.tm_mday      << " ";
    ss << std::setfill('0') << std::setw(2) << tm_utc_start.tm_hour      << ":";
    ss << std::setfill('0') << std::setw(2) << tm_utc_start.tm_min       << ":";
    ss << std::setfill('0') << std::setw(2) << tm_utc_start.tm_sec;

    return ss.str();
}

template<typename TF>
Interpolation_factors<TF> Timeloop<TF>::get_interpolation_factors(const std::vector<double>& timevec)
{
    // 1. Get the indexes and factors for the interpolation in time
    std::vector<unsigned long> itimevec(timevec.size());
    for (size_t t=0; t<timevec.size(); ++t)
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
        std::string msg = " Interpolation time is out of range; t0 = " + std::to_string(timevec[0]) + "; current time is " + std::to_string(time);
        throw std::runtime_error(msg);
    }
    else if (ifac.index1 == timevec.size())
    {
        std::string msg = " Interpolation time is out of range; t1 = " + std::to_string(timevec[timevec.size()-1]) + "; current time is " + std::to_string(time);
        throw std::runtime_error(msg);
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
