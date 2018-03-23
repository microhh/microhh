/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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

#include <vector>
#include <string>
#include <iostream>

#include "radiation.h"
#include "input.h"
#include "data_block.h"

namespace
{
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
}

template<typename TF>
Radiation<TF>::Radiation(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    // Read the switches from the input
    std::string swradiation_in = inputin.get_item<std::string>("radiation", "swradiation", "", "0");

    if (swradiation_in == "0")
        swradiation = Radiation_type::Disabled;
    else if (swradiation_in == "1")
        swradiation = Radiation_type::Enabled;
    else
        throw std::runtime_error("Invalid option for \"swradiation\"");

    double cp = 1004.;
    c_rrtmg_lw_init(&cp);
    ncol = 1;
    nlay = 60;
    nbndlw = 16;
}

template<typename TF>
Radiation<TF>::~Radiation()
{
}

template<typename TF>
void Radiation<TF>::init()
{
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
}

template<typename TF>
void Radiation<TF>::create()
{
    std::string block_name = "radiation.prof";
    Data_block data_block(master, block_name);
    data_block.get_vector(play, "pavel", nlay, 0, 0);
    data_block.get_vector(plev, "pz", nlay+1, 0, 0);
    data_block.get_vector(tlay, "tavel", nlay, 0, 0);
    data_block.get_vector(tlev, "tz", nlay+1, 0, 0);
    data_block.get_vector(h2ovmr, "wkl1", nlay, 0, 0);
    data_block.get_vector(co2vmr, "wkl2", nlay, 0, 0);
    data_block.get_vector(o3vmr, "wkl3", nlay, 0, 0);
    data_block.get_vector(n2ovmr, "wkl4", nlay, 0, 0);
    data_block.get_vector(ch4vmr, "wkl6", nlay, 0, 0);
    data_block.get_vector(o2vmr, "wkl7", nlay, 0, 0);

    tsfc[0] = 300;
}

template<typename TF>
void Radiation<TF>::exec()
{
    // For now...
    inflglw = 0;
    iceflglw = 0;
    liqflglw = 0;

    std::fill(emis.begin(), emis.end(), 1.); 
    // End

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
}

template class Radiation<double>;
template class Radiation<float>;
