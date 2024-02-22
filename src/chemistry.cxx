
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

//#include <cstdio>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <utility>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"
#include "netcdf_interface.h"
#include "chemistry.h"
#include "constants.h"
#include "timeloop.h"
#include "deposition.h"
#include "boundary.h"
#include "cross.h"

namespace
{
    // From microhh_root/kpp:
    #include "mhh_Parameters.h"
    #include "mhh_Global.h"
    #include "mhh_Sparse.h"
    #include "mhh_Integrator.c"      /* needs to be modified */
    #include "mhh_Function.c"        /* needs to be modified */
    #include "mhh_LinearAlgebra.c"
    #include "mhh_JacobianSP.c"
    #include "mhh_Jacobian.c"
    #include "mhh_Rates.c"

    double C[NSPEC];                         /* Concentration of all species */
    double * VAR = & C[0];
    double * FIX = & C[14];
    double RCONST[NREACT];                   /* Rate constants (global) */
    double RF[NREACT];                       /* Modified: These contain the reaction fluxes also modify call in Integrator */
    double Vdot[NVAR];                       /* Time derivative of variable species concentrations */
    double TIME;                             /* Current integration time */
    double RTOLS;                            /* (scalar) Relative tolerance */
    double TSTART;                           /* Integration start time */
    double TEND;                             /* Integration end time */
    double DT;                               /* Integration step */
    double ATOL[NVAR];                       /* Absolute tolerance */
    double RTOL[NVAR];                       /* Relative tolerance */
    double STEPMIN;                          /* Lower bound for integration step */
    double STEPMAX;                          /* Upper bound for integration step */
    double CFACTOR;                          /* Conversion factor for concentration units */


    std::pair<std::string, int> check_for_unique_time_dim(const std::map<std::string, int>& dims)
    {
        // Check for the existence of a unique time dimension.
        bool only_one_time_dim = false;
        std::string time_dim;
        int time_dim_length = 0;

        for (auto i : dims)
        {
            if (i.first.substr(0, 9) == "time_chem")
            {
                if (only_one_time_dim)
                    throw std::runtime_error("More than one time dimensions in input");
                else
                {
                    only_one_time_dim = true;
                    time_dim = i.first;
                    time_dim_length = i.second;
                }
            }
        }

        return std::make_pair(time_dim, time_dim_length);
    }



    template<typename TF>
    void pss(
            TF* restrict thno3, const TF* const restrict hno3,
            TF* restrict tn2o5, const TF* const restrict n2o5,
            TF* restrict th2o2, const TF* const restrict h2o2,
            TF* restrict tco, const TF* const restrict co,
            TF* restrict thcho, const TF* const restrict hcho,
            TF* restrict trooh, const TF* const restrict rooh,
            TF* restrict tc3h6, const TF* const restrict c3h6,
            TF* restrict tro2, const TF* const restrict ro2,
            TF* restrict tho2, const TF* const restrict ho2,
            TF* restrict to3, const TF* const restrict o3,
            TF* restrict tno, const TF* const restrict no,
            TF* restrict tno3, const TF* const restrict no3,
            TF* restrict tno2, const TF* const restrict no2,
            TF* restrict toh, const TF* const restrict oh,
            const TF* const restrict jval,
            const TF* const restrict emval,
            const TF* const restrict vdo3,
            const TF* const restrict vdno,
            const TF* const restrict vdno2,
            const TF* const restrict vdhno3,
            const TF* const restrict vdh2o2,
            const TF* const restrict vdrooh,
            const TF* const restrict vdhcho,
            const TF* restrict qt,
            const TF* restrict temp,
            const TF* const restrict dzi,
            const TF* const restrict rhoref,
            TF* restrict rfa,
            TF& trfa,
            const TF dt,
            const TF sdt,
            const TF switch_dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int Pj_o31d = 0;
        const int Pj_h2o2 = 1;
        const int Pj_no2  = 2;
        const int Pj_no3  = 3;
        const int Pj_n2o5 = 4;
        const int Pj_ch2or = 5;
        const int Pj_ch2om = 6;
        const int Pj_ch3o2h = 7;

        const TF RTOLS = 1e-3;

        const TF xmh2o = 18.015265;
        const TF xmh2o_i = TF(1) / xmh2o;
        const TF xmair = 28.9647;       // Molar mass of dry air  [kg kmol-1]
        const TF xmair_i = TF(1) / xmair;
        const TF Na = 6.02214086e23; // Avogadros number [molecules mol-1]

        TF C_M = 2.55e19;  // because of scope KPP generated routines
        TF VAR0[NVAR] ;
        TF erh = TF(0.0);
        TF eno = TF(0.0);

        for(int i=0; i<NVAR; i++)
        {
            RTOL[i] = RTOLS;
            ATOL[i] = 1.0;
            //tscale[i] = 1e20;
        }

        TF STEPMIN = 0.01;
        TF STEPMAX = 90;

        VAR = &C[0];
        FIX = &C[14];

        // Update the time integration of the reaction fluxes with the full timestep on first RK3 step
        if (abs(sdt/dt - 1./3.) < 1e-5) trfa += dt;

        for (int k=kstart; k<kend; ++k)
        {
            C_M = TF(1e-3) * rhoref[k] * Na * xmair_i;   // molecules/cm3 for chmistry!

            // From ppb (units mixing ratio) to molecules/cm3 --> changed: now mol/mol unit for transported tracers:
            const TF CFACTOR = C_M;
            const TF C_H2 = TF(500e-9) * CFACTOR ; // 500 ppb --> #/cm3
            const TF sdt_cfac_i = TF(1) / (sdt * CFACTOR);

            if (k==kstart)
            {
                // Emission/deposition fluxes:
                //erh      = (TF)0.1*CFACTOR/dz[k];
                //eno      = (TF)0.1*CFACTOR/dz[k];
                erh = TF(0) * CFACTOR * dzi[k];
                eno = TF(0) * CFACTOR * dzi[k];
            }
            else
            {
                erh = TF(0);
                eno = TF(0);
            }

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    const int ij = i + j*jstride;

                    // kg/kg --> molH2O/molAir --*C_M--> molecules/cm3 limit to 1 molecule/cm3 to avoid error usr_HO2_HO2
                    const TF C_H2O = std::max(qt[ijk] * xmair * C_M * xmh2o_i, TF(1));
                    const TF TEMP = temp[ijk];

                    // Convert to molecules per cm3 and add tendenccies of other processes.
                    VAR[ind_HNO3] = std::max((hno3[ijk] + thno3[ijk] * sdt) * CFACTOR, TF(0));
                    VAR[ind_N2O5] = std::max((n2o5[ijk] + tn2o5[ijk] * sdt) * CFACTOR, TF(0));
                    VAR[ind_H2O2] = std::max((h2o2[ijk] + th2o2[ijk] * sdt) * CFACTOR, TF(0));
                    VAR[ind_CO  ] = std::max((co[ijk]   + tco[ijk]   * sdt) * CFACTOR, TF(0));
                    VAR[ind_HCHO] = std::max((hcho[ijk] + thcho[ijk] * sdt) * CFACTOR, TF(0));
                    VAR[ind_ROOH] = std::max((rooh[ijk] + trooh[ijk] * sdt) * CFACTOR, TF(0));
                    VAR[ind_RH]   = std::max((c3h6[ijk] + tc3h6[ijk] * sdt) * CFACTOR, TF(0));
                    VAR[ind_RO2]  = std::max((ro2[ijk]  + tro2[ijk]  * sdt) * CFACTOR, TF(0));
                    VAR[ind_HO2]  = std::max((ho2[ijk]  + tho2[ijk]  * sdt) * CFACTOR, TF(0));
                    VAR[ind_O3]   = std::max((o3[ijk]   + to3[ijk]   * sdt) * CFACTOR, TF(0));
                    VAR[ind_NO]   = std::max((no[ijk]   + tno[ijk]   * sdt) * CFACTOR, TF(0));
                    VAR[ind_NO3]  = std::max((no3[ijk]  + tno3[ijk]  * sdt) * CFACTOR, TF(0));
                    VAR[ind_NO2]  = std::max((no2[ijk]  + tno2[ijk]  * sdt) * CFACTOR, TF(0));
                    VAR[ind_OH]   = std::max((oh[ijk]   + toh[ijk]   * sdt) * CFACTOR, TF(0));

                    RCONST[0 ] = ARR3(TF(1.7E-12), TF(-940), TEMP);
                    RCONST[1 ] = ARR3(TF(1.E-14), TF(-490), TEMP);
                    RCONST[2 ] = ARR3(TF(4.8E-11), TF(250), TEMP);
                    RCONST[3 ] = EPR(TF(3.E-13), TF(460), TF(2.1E-33), TF(920), TF(1.4E-21), TF(2200), C_M, C_H2O, TEMP);
                    RCONST[4 ] = ARR3(TF(2.9E-12), TF(-160), TEMP);
                    RCONST[5 ] = ARR3(TF(2.8E-12), TF(-1800), TEMP) * C_H2;
                    RCONST[6 ] = ARR3(TF(3.E-12), TF(-1500), TEMP);
                    RCONST[7 ] = ARR3(TF(1.4E-13), TF(-2470), TEMP);
                    RCONST[8 ] = ARR3(TF(1.8E-11), TF(110), TEMP);
                    RCONST[9 ] = TROE_ifs(TF(3.6E-30), TF(4.1), TF(1.9E-12), TF(-0.2), TF(10), C_M, TEMP);
                    RCONST[10] = TROE_ifs2(TF(1.3E-3), TF(-3.5), TF(9.7E14), TF(0.1), TF(10), C_M, TF(-11000), TF(-11080), TEMP);
                    RCONST[11] = ARR3(TF(3.3E-12), TF(270), TEMP);
                    RCONST[12] = TROE_no2oh(TF(3.2E-30), TF(4.5), TF(3.E-11), TF(10), C_M, TEMP);
                    RCONST[13] = TF(4e-12);
                    RCONST[14] = RK28(TF(2.4E-14), TF(460), TF(6.51E-34), TF(1335), TF(2.69E-17), TF(2199), C_M, TEMP);
                    RCONST[15] = TF(0.0004);
                    RCONST[16] = ARR3(TF(2.45E-12), TF(-1775), TEMP);
                    RCONST[17] = ARR3(TF(3.8E-13), TF(780), TEMP) * (TF(1)-(TF(1)/(TF(1) + ARR3(TF(498), TF(-1160), TEMP))));
                    RCONST[18] = ARR3(TF(3.8E-13), TF(780), TEMP) * (TF(1)/(TF(1.) + ARR3(TF(498.), TF(-1160), TEMP)));
                    RCONST[19] = ARR3(TF(2.8E-12), TF(300), TEMP);
                    RCONST[20] = TF(1.2e-12);
                    RCONST[21] = ARR3(TF(3.8E-12), TF(200), TEMP);
                    RCONST[22] = ARR3(TF(5.5E-12), TF(125), TEMP);
                    RCONST[23] = TF(5.8e-16);
                    RCONST[24] = TROE_cooh(TF(5.9E-33), TF(1.4), TF(1.1E-12), TF(-1.3), TF(1.5E-13), TF(-0.6), TF(2.1E9), TF(-6.1), TF(0.6), C_M, TEMP);
                    RCONST[25] = ARR3(TF(9.5E-14), TF(390), TEMP);
                    RCONST[26] = ARR3(TF(5.5E-15), TF(-1880), TEMP);
                    RCONST[27] = k3rd_iupac(TF(8.6E-27), TF(3.5), TF(3.E-11), TF(1), TF(0.6), C_M, TF(0.5), TEMP);
                    RCONST[28] = ARR3(TF(4.6E-13), TF(-1155), TEMP);
                    RCONST[29] = usr_O3_hv_H2O(TEMP, C_M, C_H2O, jval[Pj_o31d]);
                    RCONST[30] = jval[Pj_no2];
                    RCONST[31] = jval[Pj_n2o5];
                    RCONST[32] = jval[Pj_no3];
                    RCONST[33] = jval[Pj_ch3o2h];
                    RCONST[34] = jval[Pj_ch2om];
                    RCONST[35] = jval[Pj_ch2or];
                    RCONST[36] = jval[Pj_h2o2];
                    RCONST[37] = eno;
                    RCONST[38] = erh;

                    if (k==kstart)
                    {
                        RCONST[39] = vdo3[ij]   * dzi[k];
                        RCONST[40] = vdno[ij]   * dzi[k];
                        RCONST[41] = vdno2[ij]  * dzi[k];
                        RCONST[42] = vdhno3[ij] * dzi[k];
                        RCONST[43] = vdh2o2[ij] * dzi[k];
                        RCONST[44] = vdhcho[ij] * dzi[k];
                        RCONST[45] = vdrooh[ij] * dzi[k];
                    }
                    else
                    {
                        RCONST[39] = TF(0.0);
                        RCONST[40] = TF(0.0);
                        RCONST[41] = TF(0.0);
                        RCONST[42] = TF(0.0);
                        RCONST[43] = TF(0.0);
                        RCONST[44] = TF(0.0);
                        RCONST[45] = TF(0.0);
                    }

                    FIX[0] = TF(1800e-9) * CFACTOR;  // methane concentation
                    FIX[1] = C_M;                    // air density
                    FIX[2] = TF(1);                  // species added to emit

                    Fun(VAR, FIX, RCONST, Vdot, RF);

                    // Get statistics for reaction fluxes:
                    if (abs(sdt/dt - 1./3.) < 1e-5)
                    {
                        for (int l=0; l<NREACT; ++l)
                            rfa[(k-kstart)*NREACT+l] +=  RF[l]*dt;    // take the first evaluation in the RK3 steps, but with full time step.
                    }

                    WCOPY(NVAR, VAR, 1, VAR0, 1);
                    INTEGRATE(TF(0), sdt);  //brings VAR0 --> VAR, with timestep sdt

                    thno3[ijk] += (VAR[ind_HNO3] - VAR0[ind_HNO3]) * sdt_cfac_i;
                    tn2o5[ijk] += (VAR[ind_N2O5] - VAR0[ind_N2O5]) * sdt_cfac_i;
                    th2o2[ijk] += (VAR[ind_H2O2] - VAR0[ind_H2O2]) * sdt_cfac_i;
                    tco[ijk] +=   (VAR[ind_CO  ] - VAR0[ind_CO  ]) * sdt_cfac_i;
                    thcho[ijk] += (VAR[ind_HCHO] - VAR0[ind_HCHO]) * sdt_cfac_i;
                    trooh[ijk] += (VAR[ind_ROOH] - VAR0[ind_ROOH]) * sdt_cfac_i;
                    tc3h6[ijk] += (VAR[ind_RH  ] - VAR0[ind_RH  ]) * sdt_cfac_i;
                    tro2[ijk] +=  (VAR[ind_RO2 ] - VAR0[ind_RO2 ]) * sdt_cfac_i;
                    tho2[ijk] +=  (VAR[ind_HO2 ] - VAR0[ind_HO2 ]) * sdt_cfac_i;
                    to3[ijk] +=   (VAR[ind_O3  ] - VAR0[ind_O3  ]) * sdt_cfac_i;
                    tno[ijk] +=   (VAR[ind_NO  ] - VAR0[ind_NO  ]) * sdt_cfac_i;
                    tno3[ijk] +=  (VAR[ind_NO3 ] - VAR0[ind_NO3 ]) * sdt_cfac_i;
                    tno2[ijk] +=  (VAR[ind_NO2 ] - VAR0[ind_NO2 ]) * sdt_cfac_i;
                    toh[ijk] +=   (VAR[ind_OH  ] - VAR0[ind_OH  ]) * sdt_cfac_i;

                    // tscale[0] = std::min(tscale[0],ABS(h2o2[ijk])/ABS(th2o2[ijk]));
                    // tscale[1] = std::min(tscale[1],ABS(ch4[ijk])/ABS(tch4[ijk]));
                    // tscale[2] = std::min(tscale[2],ABS(n2o5[ijk])/ABS(tn2o5[ijk]));
                    // tscale[3] = std::min(tscale[3],ABS(hald[ijk])/ABS(thald[ijk]));
                    // tscale[4] = std::min(tscale[4],ABS(co[ijk])/ABS(tco[ijk]));
                    // tscale[5] = std::min(tscale[5],ABS(hcho[ijk])/ABS(thcho[ijk]));
                    // tscale[6] = std::min(tscale[6],ABS(isopooh[ijk])/ABS(tisopooh[ijk]));
                    // tscale[7] = std::min(tscale[7],ABS(isop[ijk])/ABS(tisop[ijk]));
                    // tscale[8] = std::min(tscale[8],ABS(mvkmacr[ijk])/ABS(tmvkmacr[ijk]));
                    // tscale[9] = std::min(tscale[9],ABS(xo2[ijk])/ABS(txo2[ijk]));
                    // tscale[10] = std::min(tscale[10],ABS(isopao2[ijk])/ABS(tisopao2[ijk]));
                    // tscale[11] = std::min(tscale[11],ABS(no2[ijk])/ABS(tno2[ijk]));
                    // tscale[12] = std::min(tscale[12],ABS(o3[ijk])/ABS(to3[ijk]));
                    // tscale[13] = std::min(tscale[13],ABS(no[ijk])/ABS(tno[ijk]));
                    // tscale[14] = std::min(tscale[14],ABS(ch3o2[ijk])/ABS(tch3o2[ijk]));
                    // tscale[15] = std::min(tscale[15],ABS(isopbo2[ijk])/ABS(tisopbo2[ijk]));
                    // tscale[16] = std::min(tscale[16],ABS(no3[ijk])/ABS(tno3[ijk]));
                    // tscale[17] = std::min(tscale[17],ABS(ho2[ijk])/ABS(tho2[ijk]));
                    // tscale[18] = std::min(tscale[18],ABS(oh[ijk])/ABS(toh[ijk]));
                } // i
        } // k
    }
}

template<typename TF>
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();

    sw_chemistry = inputin.get_item<bool>("chemistry", "swchemistry", "", false);

    if (!sw_chemistry)
        return;

    deposition = std::make_shared<Deposition <TF>>(masterin, gridin, fieldsin, inputin);
}

template <typename TF>
Chemistry<TF>::~Chemistry()
{
}

template<typename TF>
void Chemistry<TF>::exec_stats(const int iteration, const double time, Stats<TF>& stats)
{
    if (!sw_chemistry or stats.get_switch())
        return;

    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    auto& gd = grid.get_grid_data();

    if (iteration != 0)   // this does not make sense for first step = t=0.
    {
        // add deposition velocities to statistics:
        stats.calc_stats_2d("vdo3"   , vdo3,   no_offset);
        stats.calc_stats_2d("vdno"   , vdno,   no_offset);
        stats.calc_stats_2d("vdno2"  , vdno2,  no_offset);
        stats.calc_stats_2d("vdhno3" , vdhno3, no_offset);
        stats.calc_stats_2d("vdh2o2" , vdh2o2, no_offset);
        stats.calc_stats_2d("vdrooh" , vdrooh, no_offset);
        stats.calc_stats_2d("vdhcho" , vdhcho, no_offset);

        // sum of all PEs:
        // printf("trfa: %13.4e iteration: %i time: %13.4e \n", trfa,iteration,time);
        master.sum(rfa.data(),NREACT*gd.ktot);
        for (int l=0; l<NREACT*gd.ktot; ++l)
            rfa[l] /= (trfa*gd.itot*gd.jtot);    // mean over the horizontal plane in molecules/(cm3 * s)

        // Put the data into the NetCDF file.
        const std::vector<int> time_index{statistics_counter};

        // Write the time and iteration number.
        m.time_var->insert(time     , time_index);
        m.iter_var->insert(iteration, time_index);

        const std::vector<int> time_rfaz_index = {statistics_counter, 0};

        m.profs.at("chem_budget").data = rfa;

        const int ksize = NREACT*gd.ktot;
        std::vector<int> time_rfaz_size  = {1, ksize};
        std::vector<TF> prof_nogc(
            m.profs.at("chem_budget").data.begin() ,
            m.profs.at("chem_budget").data.begin() + ksize);

        m.profs.at("chem_budget").ncvar.insert(prof_nogc, time_rfaz_index, time_rfaz_size);

        // Synchronize the NetCDF file.
        m.data_file->sync();
        // Increment the statistics index.
        ++statistics_counter;

    }

    // (re-)intialize statistics
    for (int l=0; l<NREACT*gd.ktot; ++l)
        rfa[l] = 0.0;
    trfa = (TF) 0.0;
}

template <typename TF>
void Chemistry<TF>::init(Input& inputin)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    switch_dt = inputin.get_item<TF>("chemistry", "switch_dt", "", (TF)1e5);
    statistics_counter = 0;

    // initialize 2D deposition arrays:
    vdo3.resize(gd.ijcells);
    vdno.resize(gd.ijcells);
    vdno2.resize(gd.ijcells);
    vdhno3.resize(gd.ijcells);
    vdh2o2.resize(gd.ijcells);
    vdrooh.resize(gd.ijcells);
    vdhcho.resize(gd.ijcells);

    // initialize deposition routine:
    deposition-> init(inputin);

    // fill deposition with standard values:
    std::fill(vdo3.begin(), vdo3.end(), deposition-> get_vd("o3"));
    std::fill(vdno.begin(), vdno.end(), deposition-> get_vd("no"));
    std::fill(vdno2.begin(), vdno2.end(), deposition-> get_vd("no2"));
    std::fill(vdhno3.begin(), vdhno3.end(), deposition-> get_vd("hno3"));
    std::fill(vdh2o2.begin(), vdh2o2.end(), deposition-> get_vd("h2o2"));
    std::fill(vdrooh.begin(), vdrooh.end(), deposition-> get_vd("rooh"));
    std::fill(vdhcho.begin(), vdhcho.end(), deposition-> get_vd("hcho"));

    master.print_message("Deposition arrays initialized, e.g. with vdo3 = %13.5e m/s \n", deposition-> get_vd("o3"));
}

template <typename TF>
void Chemistry<TF>::create(
        const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();
    int iotime = timeloop.get_iotime();

    Netcdf_group& group_nc = input_nc.get_group("timedep_chem");
    int time_dim_length;
    std::string time_dim;

    for (std::string varname : jname)    // check dimensions:
    {
        std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
        std::pair<std::string, int> unique_time = check_for_unique_time_dim(dims);
        time_dim = unique_time.first;
        time_dim_length = unique_time.second;
        time.resize(time_dim_length);
    }

    for (std::string varname : ename)    // check dimension also of emissions
    {
        std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
        std::pair<std::string, int> unique_time = check_for_unique_time_dim(dims);
        time_dim = unique_time.first;
        time_dim_length = unique_time.second;
        time.resize(time_dim_length);
    }

    jo31d.resize(time_dim_length);
    jh2o2.resize(time_dim_length);
    jno2.resize(time_dim_length);
    jno3.resize(time_dim_length);
    jn2o5.resize(time_dim_length);
    jch2or.resize(time_dim_length);
    jch2om.resize(time_dim_length);
    jch3o2h.resize(time_dim_length);
    emi_isop.resize(time_dim_length);
    emi_no.resize(time_dim_length);

    group_nc.get_variable(time, time_dim, {0}, {time_dim_length});
    group_nc.get_variable(jo31d, jname[0],  {0}, {time_dim_length});
    group_nc.get_variable(jh2o2, jname[1],  {0}, {time_dim_length});
    group_nc.get_variable(jno2, jname[2],  {0}, {time_dim_length});
    group_nc.get_variable(jno3, jname[3],  {0}, {time_dim_length});
    group_nc.get_variable(jn2o5, jname[4],  {0}, {time_dim_length});
    group_nc.get_variable(jch2or, jname[5],  {0}, {time_dim_length});
    group_nc.get_variable(jch2om, jname[6],  {0}, {time_dim_length});
    group_nc.get_variable(jch3o2h, jname[7],  {0}, {time_dim_length});
    group_nc.get_variable(emi_isop, ename[0],  {0}, {time_dim_length});
    group_nc.get_variable(emi_no,   ename[1],  {0}, {time_dim_length});

    // Store output of averaging.
    rfa.resize(NREACT*gd.ktot);
    for (int l=0;l<NREACT*gd.ktot;++l)
        rfa[l] = 0.0;
    trfa = (TF)0.0;

    if (stats.get_switch())
    {
        // Stats:
        const std::string group_name = "default";
        const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
        const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
        const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};

        std::stringstream filename;
        filename << sim_name << "." << "chemistry" << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";

        // Create new NetCDF file in Mask<TF> m
        m.data_file = std::make_unique<Netcdf_file>(master, filename.str(), Netcdf_mode::Create);

        // Create dimensions.
        m.data_file->add_dimension("z", gd.kmax);
        m.data_file->add_dimension("zh", gd.kmax+1);
        m.data_file->add_dimension("rfaz", NREACT*gd.ktot);
        m.data_file->add_dimension("ijcells",gd.ijcells);
        m.data_file->add_dimension("time");

        // Create variables belonging to dimensions.
        Netcdf_handle& iter_handle =
                m.data_file->group_exists("default") ? m.data_file->get_group("default") : m.data_file->add_group("default");

        m.iter_var = std::make_unique<Netcdf_variable<int>>(iter_handle.add_variable<int>("iter", {"time"}));
        m.iter_var->add_attribute("units", "-");
        m.iter_var->add_attribute("long_name", "Iteration number");

        m.time_var = std::make_unique<Netcdf_variable<TF>>(m.data_file->template add_variable<TF>("time", {"time"}));
        if (timeloop.has_utc_time())
            m.time_var->add_attribute("units", "seconds since " + timeloop.get_datetime_utc_start_string());
        else
            m.time_var->add_attribute("units", "seconds since start");
        m.time_var->add_attribute("long_name", "Time");

        Netcdf_variable<TF> z_var = m.data_file->template add_variable<TF>("z", {"z"});
        z_var.add_attribute("units", "m");
        z_var.add_attribute("long_name", "Full level height");

        Netcdf_variable<TF> zh_var = m.data_file->template add_variable<TF>("zh", {"zh"});
        zh_var.add_attribute("units", "m");
        zh_var.add_attribute("long_name", "Half level height");

        std::string name = "chem_budget";
        std::string longname = "chemistry budget per layer";
        std::string unit = "molecules cm-3 s-1";
        Netcdf_variable<TF> rfaz_var = m.data_file->template add_variable<TF>("rfaz", {"rfaz"});
        rfaz_var.add_attribute("units", unit);
        rfaz_var.add_attribute("long_name", longname);

        // add a profile of reaction rates x z
        Level_type level =  Level_type::Full;

        Netcdf_handle& handle =
                m.data_file->group_exists("default") ? m.data_file->get_group("default") : m.data_file->add_group("default");
        Prof_var<TF> tmp{handle.add_variable<TF>(name, {"time", "rfaz"}), std::vector<TF>(gd.ktot*NREACT), level};
        m.profs.emplace(
                std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(tmp)));

        m.profs.at(name).ncvar.add_attribute("units", unit);
        m.profs.at(name).ncvar.add_attribute("long_name", longname);

        // Save the grid variables.
        std::vector<TF> z_nogc (gd.z. begin() + gd.kstart, gd.z. begin() + gd.kend  );
        std::vector<TF> zh_nogc(gd.zh.begin() + gd.kstart, gd.zh.begin() + gd.kend+1);
        z_var .insert( z_nogc, {0});
        zh_var.insert(zh_nogc, {0});

        // Synchronize the NetCDF file.
        m.data_file->sync();

        m.nmask. resize(gd.kcells);
        m.nmaskh.resize(gd.kcells);

        // add the deposition-velocity timeseries in deposition group statistics
        const std::string group_named = "deposition";

        // used in chemistry:
        stats.add_time_series("vdo3", "O3 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdno", "NO deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdno2", "NO2 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdhno3", "HNO3 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdh2o2", "H2O2 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdrooh", "ROOH deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdhcho", "HCHO deposition velocity", "m s-1", group_named);
    }

    // add cross-sections
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {"vdo3", "vdno", "vdno2", "vdhno3", "vdh2o2", "vdrooh", "vdhcho"};
        cross_list = cross.get_enabled_variables(allowed_crossvars);

        // `deposition->create()` only creates cross-sections.
        deposition->create(stats, cross);
    }
}

template<typename TF>
void Chemistry<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    for (auto& it : cross_list)
    {
        if (it == "vdo3")
            cross.cross_plane(vdo3.data(), "vdo3", iotime);
        else if (it == "vdno")
            cross.cross_plane(vdno.data(), "vdno", iotime);
        else if (it == "vdno2")
            cross.cross_plane(vdno2.data(), "vdno2", iotime);
        else if (it == "vdhno3")
            cross.cross_plane(vdhno3.data(), "vdhno3", iotime);
        else if (it == "vdh2o2")
            cross.cross_plane(vdh2o2.data(), "vdh2o2", iotime);
        else if (it == "vdrooh")
            cross.cross_plane(vdrooh.data(), "vdrooh", iotime);
        else if (it == "vdhcho")
            cross.cross_plane(vdhcho.data(), "vdhcho", iotime);
    }

    // see if to write per tile:
    deposition->exec_cross(cross, iotime);
}

template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary)
{
    if (!sw_chemistry)
        return;

    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);

    jval[0] = ifac.fac0 * jo31d[ifac.index0] + ifac.fac1 * jo31d[ifac.index1];
    jval[1] = ifac.fac0 * jh2o2[ifac.index0] + ifac.fac1 * jh2o2[ifac.index1];
    jval[2] = ifac.fac0 * jno2[ifac.index0]  + ifac.fac1 * jno2[ifac.index1];
    jval[3] = ifac.fac0 * jno3[ifac.index0] + ifac.fac1 * jno3[ifac.index1];
    jval[4] = ifac.fac0 * jn2o5[ifac.index0] + ifac.fac1 * jn2o5[ifac.index1];
    jval[5] = ifac.fac0 * jch2or[ifac.index0] + ifac.fac1 * jch2or[ifac.index1];
    jval[6] = ifac.fac0 * jch2om[ifac.index0] + ifac.fac1 * jch2om[ifac.index1];
    jval[7] = ifac.fac0 * jch3o2h[ifac.index0] + ifac.fac1 * jch3o2h[ifac.index1];
    emval[0] = ifac.fac0 * emi_isop[ifac.index0] + ifac.fac1 * emi_isop[ifac.index1];
    emval[1] = ifac.fac0 * emi_no[ifac.index0] + ifac.fac1 * emi_no[ifac.index1];

    deposition->update_time_dependent(
            timeloop,
            boundary,
            vdo3.data(),
            vdno.data(),
            vdno2.data(),
            vdhno3.data(),
            vdh2o2.data(),
            vdrooh.data(),
            vdhcho.data());
}


#ifndef USECUDA
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo,double sdt,double dt)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    auto tmp = fields.get_tmp();
    thermo.get_thermo_field(*tmp, "T", true, false);

    pss<TF>(
        fields.st.at("hno3")->fld.data(), fields.sp.at("hno3")->fld.data(),
        fields.st.at("n2o5")->fld.data(), fields.sp.at("n2o5")->fld.data(),
        fields.st.at("h2o2")->fld.data(), fields.sp.at("h2o2")->fld.data(),
        fields.st.at("co")  ->fld.data(), fields.sp.at("co")->fld.data(),
        fields.st.at("hcho")->fld.data(), fields.sp.at("hcho")->fld.data(),
        fields.st.at("rooh")->fld.data(), fields.sp.at("rooh")->fld.data(),
        fields.st.at("c3h6")->fld.data(), fields.sp.at("c3h6")->fld.data(),
        fields.st.at("ro2") ->fld.data(), fields.sp.at("ro2")->fld.data(),
        fields.st.at("ho2") ->fld.data(), fields.sp.at("ho2")->fld.data(),
        fields.st.at("o3")  ->fld.data(), fields.sp.at("o3")->fld.data(),
        fields.st.at("no")  ->fld.data(), fields.sp.at("no")->fld.data(),
        fields.st.at("no3") ->fld.data(), fields.sp.at("no3")->fld.data(),
        fields.st.at("no2") ->fld.data(), fields.sp.at("no2")->fld.data(),
        fields.st.at("oh")  ->fld.data(), fields.sp.at("oh")->fld.data(),
        jval, emval,
        vdo3.data(),
        vdno.data(),
        vdno2.data(),
        vdhno3.data(),
        vdh2o2.data(),
        vdrooh.data(),
        vdhcho.data(),
        fields.sp.at("qt")->fld.data(),
        tmp->fld.data(),
        gd.dzi.data(),
        fields.rhoref.data(),
        rfa.data(),
        trfa,
        dt, sdt, switch_dt,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);

    fields.release_tmp(tmp);

    //isop_stat<TF>(
    //      fields.st.at("isop")->fld.data(), fields.sp.at("isop")->fld.data(),
    //      gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
    //      gd.icells, gd.ijcells);
}
#endif

template class Chemistry<double>;
//:template class Chemistry<float>;
