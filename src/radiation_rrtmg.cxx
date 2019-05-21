/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include "radiation.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "radiation_rrtmg.h"

namespace
{
    /*
    // Wrapper functions to the RRTMG long wave kernels.
    extern "C"
    {
        void c_rrtmg_lw_init(double *cpdair);
        void c_rrtmg_lw (
                int *ncol    ,int *nlay    ,int *icld    ,int *idrv    ,
                double *play    ,double *plev    ,double *tlay    ,double *tlev    ,double *tsfc    ,
                double *h2ovmr  ,double *o3vmr   ,double *co2vmr  ,double *ch4vmr  ,double *n2ovmr  ,double *o2vmr,
                double *cfc11vmr,double *cfc12vmr,double *cfc22vmr,double *ccl4vmr ,double *emis    ,
                int *inflglw ,int *iceflglw,int *liqflglw,double *cldfr   ,
                double *taucld  ,double *cicewp  ,double *cliqwp  ,double *reice   ,double *reliq   ,
                double *tauaer  ,
                double *uflx    ,double *dflx    ,double *hr      ,double *uflxc   ,double *dflxc,  double *hrc,
                double *duflx_dt,double *duflxc_dt );
    }
    */
}

template<typename TF>
Radiation_rrtmg<TF>::Radiation_rrtmg(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Radiation<TF>(masterin, gridin, fieldsin, inputin)
{
    swradiation = "1";

    // double cp = 1004.;
    // c_rrtmg_lw_init(&cp);
    // ncol = 1;
    // nlay = 60;
    // nbndlw = 16;
}

template<typename TF>
Radiation_rrtmg<TF>::~Radiation_rrtmg()
{
}

template<typename TF>
void Radiation_rrtmg<TF>::init()
{
    /*
    play.resize(ncol*nlay);     // (ncol, nlay)
    plev.resize(ncol*(nlay+1)); // (ncol, nlay+1)
    tlay.resize(ncol*nlay);     // (ncol, nlay)
    tlev.resize(ncol*(nlay+1)); // (ncol, nlay+1)

    tsfc.resize(ncol); // (ncol)

    h2ovmr.resize(ncol*nlay); // (ncol, nlay)
    o3vmr.resize(ncol*nlay);  // (ncol, nlay)
    co2vmr.resize(ncol*nlay); // (ncol, nlay)
    ch4vmr.resize(ncol*nlay); // (ncol, nlay)
    n2ovmr.resize(ncol*nlay); // (ncol, nlay)
    o2vmr.resize(ncol*nlay);  // (ncol, nlay)

    cfc11vmr.resize(ncol*nlay); // (ncol, nlay)
    cfc12vmr.resize(ncol*nlay); // (ncol, nlay)
    cfc22vmr.resize(ncol*nlay); // (ncol, nlay)
    ccl4vmr.resize(ncol*nlay);  // (ncol, nlay)
    emis.resize(ncol*nbndlw);   // (ncol, nbndlw)

    cldfr.resize(ncol*nlay);  // (ncol, nlay)
    cicewp.resize(ncol*nlay); // (ncol, nlay)
    cliqwp.resize(ncol*nlay); // (ncol, nlay)
    reice.resize(ncol*nlay);  // (ncol, nlay)
    reliq.resize(ncol*nlay);  // (ncol, nlay)
    taucld.resize(nbndlw*ncol*nlay); // (nbndlw, ncol, nlay)
    tauaer.resize(ncol*nlay*nbndlw); // (ncol, nlay, nbndlw)

    // OUTPUT
    uflx.resize(ncol*(nlay+1));      // (ncol, nlay+1)
    dflx.resize(ncol*(nlay+1));      // (ncol, nlay+1)
    hr.resize(ncol*nlay);            // (ncol, nlay)
    uflxc.resize(ncol*(nlay+1));     // (ncol, nlay+1)
    dflxc.resize(ncol*(nlay+1));     // (ncol, nlay+1)
    hrc.resize(ncol*nlay);           // (ncol, nlay)
    duflx_dt.resize(ncol*(nlay+1));  // (ncol, nlay+1)
    duflxc_dt.resize(ncol*(nlay+1)); // (ncol, nlay+1)
    */
}

template<typename TF>
void Radiation_rrtmg<TF>::create(Thermo<TF>& thermo,Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
    /*
    std::string block_name = "radiation.prof";
    Data_block data_block(master, block_name);

    std::vector<TF> p_rad(nlay);
    data_block.get_vector(p_rad, "pavel", nlay, 0, 0);
    std::vector<TF> ph_rad(nlay+1);
    data_block.get_vector(ph_rad, "pz", nlay+1, 0, 0);

    // Convert the data to
    for (auto& d : p_rad)
        d *= 100.;

    for (auto& d : ph_rad)
        d *= 100.;

    // Get copies of pressure profiles.
    std::vector<TF> p = thermo.get_p_vector();
    std::vector<TF> ph = thermo.get_ph_vector();

    auto& gd = grid.get_grid_data();

    // Remove the ghost cells from the vector.
    p.erase(p.begin() + gd.kend, p.end());
    ph.erase(ph.begin() + gd.kend+1, ph.end());
    p.erase(p.begin(), p.begin() + gd.kstart);
    ph.erase(ph.begin(), ph.begin() + gd.kstart);

    // Remove all elements from radiation profile vector that are in the domain.
    auto it = p_rad.begin();
    int counter = 0;
    while (*it > ph.back())
    {
        ++counter;
        ++it;
    }

    // Delete all items until first index, remove one level extra for ph
    p_rad.erase(p_rad.begin(), p_rad.begin()+counter);
    ph_rad.erase(ph_rad.begin(), ph_rad.begin()+counter+1);
    p_rad.insert(p_rad.begin(), p.begin(), p.end());
    ph_rad.insert(ph_rad.begin(), ph.begin(), ph.end());

    // CvH if the first pressure level above the domain top is too close to the half level, it could be moved...

    // Temporary...
    data_block.get_vector(play, "pavel", nlay, 0, 0);
    data_block.get_vector(plev, "pz", nlay+1, 0, 0);
    data_block.get_vector(tlay, "tavel", nlay, 0, 0);
    data_block.get_vector(tlev, "tz", nlay+1, 0, 0);

    // The specific humidity comes from climatology for now.
    data_block.get_vector(h2ovmr, "wkl1", nlay, 0, 0);

    // These profiles come from climatology
    data_block.get_vector(co2vmr, "wkl2", nlay, 0, 0);
    data_block.get_vector(o3vmr, "wkl3", nlay, 0, 0);
    data_block.get_vector(n2ovmr, "wkl4", nlay, 0, 0);
    data_block.get_vector(ch4vmr, "wkl6", nlay, 0, 0);
    data_block.get_vector(o2vmr, "wkl7", nlay, 0, 0);
    */
}

template<typename TF>
void Radiation_rrtmg<TF>::exec(Thermo<TF>& thermo, double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    /*
    // auto& gd = grid.get_grid_data();

    // For now...
    inflglw = 0;
    iceflglw = 0;
    liqflglw = 0;

    std::fill(emis.begin(), emis.end(), 1.);

    // Step 1. Get the mean absolute atmospheric temperature.
    auto T = fields.get_tmp();
    thermo.get_thermo_field(*T, "T", false, false);

    // Calculate radiative cooling only for single column.
    field3d_operators.calc_mean_profile(T->fld_mean.data(), T->fld.data());

    thermo.get_thermo_field(*T, "T_h", false, false);

    // Step 2.

    // Get absolute atmospheric and surface temperature.
    tsfc[0] = 300.;

    c_rrtmg_lw(
            &ncol    ,&nlay    ,&icld    ,&idrv    ,
            play.data()    ,plev.data()    ,tlay.data()    ,tlev.data()    ,tsfc.data()    ,
            h2ovmr.data()  ,o3vmr.data()   ,co2vmr.data()  ,ch4vmr.data()  ,n2ovmr.data()  ,o2vmr.data(),
            cfc11vmr.data(),cfc12vmr.data(),cfc22vmr.data(),ccl4vmr.data() ,emis.data()    ,
            &inflglw ,&iceflglw,&liqflglw,cldfr.data() ,
            taucld.data()  ,cicewp.data()  ,cliqwp.data()  ,reice.data()   ,reliq.data() ,
            tauaer.data()  ,
            uflx.data()    ,dflx.data()    ,hr.data()      ,uflxc.data()   ,dflxc.data(),  hrc.data(),
            duflx_dt.data(),duflxc_dt.data());

    fields.release_tmp(T);

    std::cout << "Heating rate" << std::endl;
    for (int i=0; i<nlay; ++i)
        std::cout << i       << ", "
            << play[i] << ", "
            << hr[i] << ", "
            << std::endl;

    std::cout << "Upflux/downflux rate" << std::endl;
    for (int i=0; i<nlay+1; ++i)
        std::cout << i       << ", "
            << plev[i] << ", "
            << uflx[i] << ", "
            << dflx[i] << ", "
            << std::endl;
    */
}

template<typename TF>
bool Radiation_rrtmg<TF>::check_field_exists(std::string name)
{
    return false;
}

template class Radiation_rrtmg<double>;
template class Radiation_rrtmg<float>;
