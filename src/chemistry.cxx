
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


namespace 
{
#include "../cases/jaenschwalde/include/tm5_ifs_Parameters.h"
#include "../cases/jaenschwalde/include/tm5_ifs_Global.h"
#include "../cases/jaenschwalde/include/tm5_ifs_Sparse.h"
#include "../cases/jaenschwalde/include/tm5_ifs_Integrator.c"      /* needs to be modified */
#include "../cases/jaenschwalde/include/tm5_ifs_Function.c"        /* needs to be modified */
#include "../cases/jaenschwalde/include/tm5_ifs_LinearAlgebra.c"
#include "../cases/jaenschwalde/include/tm5_ifs_JacobianSP.c"
#include "../cases/jaenschwalde/include/tm5_ifs_Jacobian.c"
	

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
        TF EPR(const TF A1, const TF C1, const TF A2, const TF C2, 
           const TF A3, const TF C3, const TF mmult, const TF ch2o,
	   const TF TEMP) 
	      {               
	      TF K1, K2, K3, EPR_p;
	      
	      K1 = (TF)A1 * exp(C1/TEMP);
	      K2 = (TF)A2 * exp(C2/TEMP) * mmult;
	      K3 = (TF)A3 * exp(C3/TEMP) * ch2o;
	      EPR_p = (K1 + K2) * (1.0 + K3);
	      
	      return (TF)EPR_p;
	      }

	template<typename TF>
	TF ARR3(const TF A0, const TF B0, const TF TEMP)
	       {
	       return A0 * exp(B0 / TEMP);
	       }


	template<typename TF>
	TF TROE_no2oh(const TF kzero, const TF mzero, const TF kinf, 
		    const TF fmulti, const TF MN2, const TF TEMP)
	       {
	       TF k0T, kinfT, znn;
	       k0T  = (kzero * pow(((TF)300./TEMP), mzero)) * MN2;
	       kinfT = kinf;
	       znn = (TF)0.75 - ((TF)1.27 * log10((TF)0.41));
	       return (k0T * kinfT) / (k0T + kinfT) * pow(fmulti, 
		      (log10((TF)0.41) / ((TF)1. + pow((log10(k0T/kinfT))/znn, 2))));
	       }

	template<typename TF>
	TF TROE_cooh(const TF kzero, const TF mzero, double kinf, 
		    const TF minf, const TF k3, const TF c3, const TF k4, 
		    const TF c4, const TF fmulti, const TF MN2, const TF TEMP)
	       {
	       TF k0T, kinfT, kval3, kval4, kcooh;
	       k0T  = (kzero * pow(((TF)300./TEMP), mzero)) * MN2;
	       kinfT = kinf * pow(((TF)300./TEMP), minf);
	       kval3 = k3 * pow(((TF)300./TEMP), c3);
	       kval4 = (k4 * pow(((TF)300./TEMP), c4)) / MN2;
	       kcooh = k0T / ((TF)1. + k0T / kinfT) * pow(fmulti, 
		     (TF)1. / ((TF)1. + pow(log10(k0T/kinfT), 2)));
	       return kcooh + (kval3 / ((TF)1. + kval3 / kval4) * pow(fmulti, 
		      (TF)1. / ((TF)1. + pow(log10(kval3/kval4), 2))));
	       }

	template<typename TF>
	TF TROE_ifs(const TF kzero, const TF mzero, const TF kinf, 
		    const TF minf, const TF fmulti, const TF MN2, const TF TEMP)
	       {
	       TF k0T, kinfT, znn_b;
	       k0T  = (kzero * pow((TEMP/(TF)300.), mzero)) * MN2;
	       kinfT = kinf * pow((TEMP/(TF)300.), minf);
	       znn_b = (TF)0.75 - ((TF)1.27 * log10((TF)0.35));
	       return (k0T * kinfT) / (k0T + kinfT) * pow(fmulti, 
		      log10((TF)0.35) / ((TF)1. + pow(log10(k0T/kinfT)/znn_b, 2)));
	       }

	template<typename TF>
	TF TROE_ifs2(const TF kzero, const TF mzero, const TF kinf, 
		    const TF minf, const TF fmulti, const TF MN2, const TF c1, 
		    const TF c2, const TF TEMP)
	       {
	       TF k0T, kinfT, znn_b;
	       k0T  = (kzero * pow((TEMP/(TF)300.), mzero)) * MN2 * exp(c1 / TEMP);
	       kinfT = kinf * pow((TEMP/(TF)300.), minf) * exp(c2 / TEMP);
	       znn_b = (TF)0.75 - ((TF)1.27 * log10((TF)0.35));
	       return (k0T * kinfT) / (k0T + kinfT) * pow(fmulti, 
		      log10((TF)0.35) / ((TF)1. + pow(log10(k0T/kinfT)/znn_b, 2)));
	       }

	template<typename TF>
	TF k3rd_iupac(TF kzero, TF mzero, TF kinf, 
		    TF minf, TF fmulti, TF MN2, TF NM, TF TEMP)
	       {
	       TF k0T;
	       TF kinfT;
	       k0T  = (kzero * pow((TEMP/(TF)300.), mzero)) * MN2;
	       kinfT = kinf * pow((TEMP/(TF)300.), minf);
	       return (k0T / ((TF)1. + (k0T/kinfT))) * 
		       pow( fmulti, (TF)1./ ((TF)1. + pow((log10(k0T/kinfT)/NM), 2) ) );
	       }

	template<typename TF>
	TF usr_O3_hv_H2O( const TF TEMP, const TF C_M, const TF C_H2O, const TF J_O1D)
	{
		TF KH2O;
		TF KN2;
		TF KO2;
		KH2O = ((TF)1.63e-10 *C_H2O * exp((TF)60.0/TEMP)  )  ;
		KN2  = ((TF)2.15e-11 * exp((TF)110.0/TEMP) *(TF)0.79*C_M) ;
		KO2  = ((TF)3.30e-11 * exp((TF)55.0 /TEMP) *(TF)0.21*C_M) ;
		return (KH2O *J_O1D) / (KH2O + KN2 + KO2);
	}           


    template<typename TF>
    void isop_stat(
            const TF* tisop, const TF* isop, 
	    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
	{
	    TF iso_mean = TF(0.0);
	    TF iso_tend = TF(0.0);
	    int cnt = 0;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
		    iso_mean += isop[ijk];
		    iso_tend += tisop[ijk];
		    cnt += 1;
		}
	    //printf("%i  %12.3e ppb  %12.3e ppb/hour \n",k,iso_mean/cnt,iso_tend*3600.0/cnt);
	}
    }

    template<typename TF>
    void pss(
            TF* restrict thno3, const TF* const restrict hno3, 
            TF* restrict tn2o5, const TF* const restrict n2o5, 
            TF* restrict th2o2, const TF* const restrict h2o2, 
            TF* restrict tco, const TF* const restrict co, 
            TF* restrict thcho, const TF* const restrict hcho, 
            TF* restrict trooh, const TF* const restrict rooh, 
            TF* restrict tc3h8, const TF* const restrict c3h8, 
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
	    TF* restrict rfa, TF& trfa,
	    const TF* restrict qt,
	    const TF* restrict Temp, const TF dt, const TF sdt, const TF switch_dt,
	    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk, const TF* const restrict dz, TF* const restrict rhoref)
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
        const TF xmair = 28.9647;       // Molar mass of dry air  [kg kmol-1]
        const TF Na    = 6.02214086e23; // Avogadros number [molecules mol-1]
	TF C_M = 2.55e19;  // because of scop KPP generated routines
	// TF tscale[NVAR] ;
	TF VAR0[NVAR] ;
        TF erh    = TF(0.0);
        TF eno    = TF(0.0);

	for( int i = 0; i < NVAR; i++ ) {
	  RTOL[i] = RTOLS;
	  ATOL[i] = 1.0;
	  // tscale[i] = 1e20;
	}
	TF STEPMIN = 0.01;
	TF STEPMAX = 90;

	VAR = &C[0];
	FIX = &C[14];
        int nkpp = 0;
	int nderiv = 0;
	// update the time integration of the reaction fluxes with the full timestep on first RK3 step
        if (abs(sdt/dt - 1./3.) < 1e-5) trfa += dt;

        for (int k=kstart; k<kend; ++k)
            {	
	    C_M = (TF)1e-3*rhoref[k]*Na/xmair;   // molecules/cm3 for chmistry!
	    const TF CFACTOR = C_M ;               // from ppb (units mixing ratio) to molecules/cm3 --> changed: now mol/mol unit for transported tracers
	    const TF C_H2 = (TF)500.0e-9*CFACTOR ;              // 500 ppb --> #/cm3
            if (k==kstart) {
		// emission/deposition fluxes:
		//erh      = (TF)0.1*CFACTOR/dz[k];
		//eno      = (TF)0.1*CFACTOR/dz[k];
		erh      = (TF)0.0*CFACTOR/dz[k];
		eno      = (TF)0.0*CFACTOR/dz[k];
	    }
            else {
		erh      = TF(0.0);
		eno      = TF(0.0);
	    }
	    TF coh = 0.0;
	    int noh = 0;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
		    const int ijk = i + j*jj + k*kk;
                    const int ij = i + j*jj;
		    const TF C_H2O = std::max(qt[ijk]*xmair*C_M/xmh2o,(TF)1.0);                   // kg/kg --> molH2O/molAir --*C_M--> molecules/cm3 limit to 1 molecule/cm3 to avoid error usr_HO2_HO2
		    const TF TEMP = Temp[ijk];
		    // convert to molecules per cm3 and add tendenccies of other processes.
                    VAR[ind_HNO3]    = std::max((hno3[ijk]  +thno3[ijk] *sdt)*CFACTOR,(TF)0.0);
		    VAR[ind_N2O5]    = std::max((n2o5[ijk]  +tn2o5[ijk] *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_H2O2]    = std::max((h2o2[ijk]  +th2o2[ijk] *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_CO  ]    = std::max((co[ijk]    +tco[ijk]   *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_HCHO]    = std::max((hcho[ijk]  +thcho[ijk] *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_ROOH]    = std::max((rooh[ijk]  +trooh[ijk] *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_RH]      = std::max((c3h8[ijk]  +tc3h8[ijk] *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_RO2]     = std::max((ro2[ijk]   +tro2[ijk]  *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_HO2]     = std::max((ho2[ijk]   +tho2[ijk]  *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_O3]      = std::max((o3[ijk]    +to3[ijk]   *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_NO]      = std::max((no[ijk]    +tno[ijk]   *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_NO3]     = std::max((no3[ijk]   +tno3[ijk]  *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_NO2]     = std::max((no2[ijk]   +tno2[ijk]  *sdt)*CFACTOR,(TF)0.0);
                    VAR[ind_OH]      = std::max((oh[ijk]    +toh[ijk]   *sdt)*CFACTOR,(TF)0.0);
		    coh += oh[ijk];
		    noh += 1;

		    RCONST[0] = (ARR3((TF)1.7E-12,(TF)-940.,TEMP));
		    RCONST[1] = (ARR3((TF)1.E-14,(TF)-490.,TEMP));
		    RCONST[2] = (ARR3((TF)4.8E-11,(TF)250.,TEMP));
		    RCONST[3] = (EPR((TF)3.E-13,(TF)460.,(TF)2.1E-33,(TF)920.,(TF)1.4E-21,(TF)2200.,C_M,C_H2O,
			     TEMP));
		    RCONST[4] = (ARR3((TF)2.9E-12,(TF)-160.,TEMP));
		    RCONST[5] = (ARR3((TF)2.8E-12,(TF)-1800.,TEMP)*C_H2);
		    RCONST[6] = (ARR3((TF)3.E-12,(TF)-1500.,TEMP));
		    RCONST[7] = (ARR3((TF)3.3E-12,(TF)270.,TEMP));
		    RCONST[8] = (TROE_no2oh((TF)3.2E-30,(TF)4.5,(TF)3.E-11,(TF)10.,C_M,TEMP));
		    RCONST[9] = (ARR3((TF)2.45E-12,(TF)-1775.,TEMP));
		    RCONST[10] = (ARR3((TF)2.8E-12,(TF)300.,TEMP));
		    RCONST[11] = (ARR3((TF)3.8E-13,(TF)780.,TEMP)*((TF)1.-((TF)1./((TF)1.+ARR3((TF)498.,(TF)-1160.,
			      TEMP)))));
		    RCONST[12] = (ARR3((TF)3.8E-13,(TF)780.,TEMP)*((TF)1./((TF)1.+ARR3((TF)498.,(TF)-1160.,
			      TEMP))));
		    RCONST[13] = (ARR3((TF)3.8E-12,(TF)200.,TEMP));
		    RCONST[14] = 2e-11;
		    RCONST[15] = (TROE_cooh((TF)5.9E-33,(TF)1.4,(TF)1.1E-12,(TF)-1.3,(TF)1.5E-13,(TF)-0.6,(TF)2.1E9,
			      (TF)-6.1,(TF)0.6,C_M,TEMP));
		    RCONST[16] = (ARR3((TF)5.5E-12,(TF)125.,TEMP));
	    	    RCONST[17] = (ARR3((TF)5.5E-15,(TF)-1880.0,TEMP));
		    RCONST[18] = (k3rd_iupac((TF)8.6E-27,(TF)3.5,(TF)3.E-11,(TF)1.,(TF)0.6,C_M,(TF)0.5,TEMP));
                    RCONST[19] = (usr_O3_hv_H2O(TEMP,C_M,C_H2O,jval[Pj_o31d]));
                    RCONST[20] = (jval[Pj_no2]);
                    RCONST[21] = (jval[Pj_ch2or]);
                    RCONST[22] = (jval[Pj_ch2om]);
		    RCONST[23] = (jval[Pj_ch3o2h]);
		    RCONST[24] = (jval[Pj_h2o2]);
		    RCONST[25] = (ARR3((TF)1.4E-13,(TF)-2470.,TEMP));
		    RCONST[26] = (ARR3((TF)1.8E-11,(TF)110.,TEMP));
		    RCONST[27] = (TROE_ifs((TF)3.6E-30,(TF)4.1,(TF)1.9E-12,(TF)-0.2,(TF)10.,C_M,TEMP));
		    RCONST[28] = ((TROE_ifs2((TF)1.3E-3,(TF)-3.5,(TF)9.7E14,(TF)0.1,(TF)10.,C_M,(TF)-11000.,
			      (TF)-11080.,TEMP)+(TF)7.348E-5));
		    RCONST[29] = (TF)4e-12;
		    RCONST[30] = (TF)1.5e-21;
		    RCONST[31] = (ARR3((TF)4.6E-13,(TF)-1155.,TEMP));
		    RCONST[32] = (jval[Pj_no3]);
		    RCONST[33] = (jval[Pj_n2o5]);
		    RCONST[34] = (erh);
		    RCONST[35] = (eno);
                    if (k==kstart) {
			    RCONST[36] = vdo3[ij]/dz[k];
			    RCONST[37] = vdno[ij]/dz[k];
			    RCONST[38] = vdno2[ij]/dz[k];
			    RCONST[39] = vdhno3[ij]/dz[k];
			    RCONST[40] = vdrooh[ij]/dz[k];
			    RCONST[41] = vdh2o2[ij]/dz[k];
			    RCONST[42] = vdhcho[ij]/dz[k];
		    }
		    else {
			    RCONST[36] = TF(0.0);
			    RCONST[37] = TF(0.0);
			    RCONST[38] = TF(0.0);
			    RCONST[39] = TF(0.0);
			    RCONST[40] = TF(0.0);
			    RCONST[41] = TF(0.0);
			    RCONST[42] = TF(0.0);
		    }
		    FIX[0] = 1800.0e-9*CFACTOR;   // methane concentation
                    FIX[1] = C_M;        // air density
		    FIX[2] = (TF)1.0;    // species added to emit   

                    // for (int ispec =0 ; ispec <14; ++ispec)
		    //    printf("%4i %13.3e ",ispec ,VAR[ispec]); 
		    //printf("\n");

		    Fun(VAR,FIX,RCONST,Vdot,RF);

		    // get statitics for reacion fluxes:
                    if (abs(sdt/dt - 1./3.) < 1e-5)
		    {
		    	for (int l=0;l<NREACT;++l) rfa[(k-kstart)*NREACT+l] +=  RF[l]*dt;    // take the first evaluation in the RK3 steps, but with full time step.
		    }

		    nkpp += 1;

		    WCOPY(NVAR,VAR,1,VAR0,1);
		    INTEGRATE(  (TF)0.0 , sdt );  //brings VAR0 --> VAR, with timestep sdt
		    thno3[ijk] +=    (VAR[ind_HNO3]-VAR0[ind_HNO3])/(sdt*CFACTOR);
		    tn2o5[ijk] +=    (VAR[ind_N2O5]-VAR0[ind_N2O5])/(sdt*CFACTOR);
		    th2o2[ijk] +=    (VAR[ind_H2O2]-VAR0[ind_H2O2])/(sdt*CFACTOR);
		    tco[ijk] +=      (VAR[ind_CO]-VAR0[ind_CO])/(sdt*CFACTOR);
		    thcho[ijk] +=    (VAR[ind_HCHO]-VAR0[ind_HCHO])/(sdt*CFACTOR);
		    trooh[ijk] +=    (VAR[ind_ROOH]-VAR0[ind_ROOH])/(sdt*CFACTOR);
		    tc3h8[ijk] +=    (VAR[ind_RH]-VAR0[ind_RH])/(sdt*CFACTOR);
		    tro2[ijk] +=     (VAR[ind_RO2]-VAR0[ind_RO2])/(sdt*CFACTOR);
		    tho2[ijk] +=     (VAR[ind_HO2]-VAR0[ind_HO2])/(sdt*CFACTOR);
		    to3[ijk] +=      (VAR[ind_O3]-VAR0[ind_O3])/(sdt*CFACTOR);
		    tno[ijk] +=      (VAR[ind_NO]-VAR0[ind_NO])/(sdt*CFACTOR);
		    tno3[ijk] +=     (VAR[ind_NO3]-VAR0[ind_NO3])/(sdt*CFACTOR);
		    tno2[ijk] +=     (VAR[ind_NO2]-VAR0[ind_NO2])/(sdt*CFACTOR);
		    toh[ijk] +=      (VAR[ind_OH]-VAR0[ind_OH])/(sdt*CFACTOR);
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
                } /* i */
	printf("%4i %13.3e %13.3e k, coh sdt \n",k,coh/noh,sdt); 
	}
    //printf("number of kpp integration %4i  number of simple derivatives %4i \n",nkpp,nderiv);
    }
}

template<typename TF>
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
	
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();
    deposition = std::make_shared<Deposition <TF>>(masterin, gridin, fieldsin, inputin);
}

template <typename TF>
Chemistry<TF>::~Chemistry()
{
}

template<typename TF>
void Chemistry<TF>::exec_stats(const int iteration, const double time, Stats<TF>& stats)
{
    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    auto& gd = grid.get_grid_data();

    for (auto& it : fields.st) if( cmap[it.first].type == Chemistry_type::disabled) return;

    /* stats.calc_stats("oh", *fields.sd.at("oh"), no_offset, no_threshold);  */

    if (iteration != 0)   // this does not make sense for first step
    {

	    // sum of all PEs: 
	    // printf("trfa: %13.4e iteration: %i time: %13.4e \n", trfa,iteration,time);
	    master.sum(rfa.data(),NREACT*gd.ktot);
	    for (int l=0;l<NREACT*gd.ktot;++l) rfa[l] /= (trfa*gd.itot*gd.jtot);    // mean over the horizontal plane in molecules/(cm3 * s)


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
    for (int l=0;l<NREACT*gd.ktot;++l) rfa[l] = 0.0;
    trfa = (TF) 0.0;
}

template <typename TF>
void Chemistry<TF>::init(Input& inputin)
{
    auto& gd = grid.get_grid_data();
    for (auto& it : fields.st)
    {
        const std::string type = inputin.get_item<std::string>("chemistry", "swchemistry", it.first, "0");
        if (type == "0")
        {
            // Cycle to avoid reading unneeded namelist options.
            continue;
        }
        else if (type == "enabled")
        {
            cmap[it.first].type = Chemistry_type::enabled;
        }
        else if (type == "disabled")
        {
            cmap[it.first].type = Chemistry_type::disabled;
        }
        else if (type == "simple")
        {
            cmap[it.first].type = Chemistry_type::simple;
        }
        else
            throw std::runtime_error("Invalid option for \"Chemistry type\"");
    }
    switch_dt = inputin.get_item<TF>("chemistry", "switch_dt","", (TF)1e5);
    statistics_counter = 0;
    // initialize 2D deposition arrays:
    vdo3.resize(gd.ijcells);
    vdno.resize(gd.ijcells);
    vdno2.resize(gd.ijcells);
    vdhno3.resize(gd.ijcells);
    vdh2o2.resize(gd.ijcells);
    vdrooh.resize(gd.ijcells);
    vdhcho.resize(gd.ijcells);
    printf("Deposition arrays created \n");
    // initialize deposition routine:
    deposition-> init(inputin);
    // fill deposition with standard values:
    std::fill(vdo3.begin(),vdo3.end(),deposition-> get_vd("o3"));
    std::fill(vdno.begin(),vdno.end(),deposition-> get_vd("no"));
    std::fill(vdno2.begin(),vdno2.end(),deposition-> get_vd("no2"));
    std::fill(vdhno3.begin(),vdhno3.end(),deposition-> get_vd("hno3"));
    std::fill(vdh2o2.begin(),vdh2o2.end(),deposition-> get_vd("h2o2"));
    std::fill(vdrooh.begin(),vdrooh.end(),deposition-> get_vd("rooh"));
    std::fill(vdhcho.begin(),vdhcho.end(),deposition-> get_vd("hcho"));
    printf("Deposition arrays initialized, e.g. with vdo3 = %13.5e m/s \n",deposition-> get_vd("o3"));




}

template <typename TF>
void Chemistry<TF>::create(const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    for (auto& it : fields.st) if( cmap[it.first].type == Chemistry_type::disabled) return;

    // nothing yet deposition->create(timeloop, sim_name, input_nc, stats);
    //

    Netcdf_group& group_nc = input_nc.get_group("timedep_chem");
    int time_dim_length;
    std::string time_dim;
    for(std::string varname : jname)    // check dimensions:
    {
	    std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
	    std::pair<std::string, int> unique_time = check_for_unique_time_dim(dims);
	    time_dim = unique_time.first;
	    time_dim_length = unique_time.second;
	    time.resize(time_dim_length);
    }
    for(std::string varname : ename)    // check dimension also of emissions
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

    // Stats:
    const std::string group_name = "default";
    const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
    const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
    const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};

    // Add the profiles to te statistics
    /* if (stats.get_switch())
    {
        stats.add_profs(*fields.sd.at("oh"), "z", stat_op_w, group_name);
    } */

//  store output of averaging
    auto& gd = grid.get_grid_data();
    rfa.resize(NREACT*gd.ktot);
    for (int l=0;l<NREACT*gd.ktot;++l) rfa[l] = 0.0;
    trfa = (TF)0.0;


    int iotime = timeloop.get_iotime();

    std::stringstream filename;
    filename << sim_name << "." << "chemistry" << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";

    // Create new NetCDF file in Mask<TF> m
    m.data_file = std::make_unique<Netcdf_file>(master, filename.str(), Netcdf_mode::Create);

    // Create dimensions.
    m.data_file->add_dimension("z",  gd.kmax);
    m.data_file->add_dimension("zh", gd.kmax+1);
    m.data_file->add_dimension("rfaz", NREACT*gd.ktot);
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
    //
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

}


template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary)
{
    for (auto& it : fields.st) if( cmap[it.first].type == Chemistry_type::disabled) return;

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
    deposition->update_time_dependent(timeloop,boundary,
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
    for (auto& it : fields.st)
        if( cmap[it.first].type == Chemistry_type::disabled) return;
    auto& gd = grid.get_grid_data();
    
    auto Temp = fields.get_tmp();
    thermo.get_thermo_field(*Temp, "T", true, false);

    pss<TF>(
	    fields.st.at("hno3")   ->fld.data(), fields.sp.at("hno3")->fld.data(), 
	    fields.st.at("n2o5")   ->fld.data(), fields.sp.at("n2o5")->fld.data(), 
	    fields.st.at("h2o2")   ->fld.data(), fields.sp.at("h2o2")->fld.data(), 
	    fields.st.at("co")     ->fld.data(), fields.sp.at("co")->fld.data(), 
	    fields.st.at("hcho")   ->fld.data(), fields.sp.at("hcho")->fld.data(), 
	    fields.st.at("rooh")   ->fld.data(), fields.sp.at("rooh")->fld.data(), 
	    fields.st.at("c3h8")   ->fld.data(), fields.sp.at("c3h8")->fld.data(), 
	    fields.st.at("ro2")    ->fld.data(), fields.sp.at("ro2")->fld.data(), 
	    fields.st.at("ho2")    ->fld.data(), fields.sp.at("ho2")->fld.data(), 
	    fields.st.at("o3")     ->fld.data(), fields.sp.at("o3")->fld.data(), 
	    fields.st.at("no")     ->fld.data(), fields.sp.at("no")->fld.data(), 
	    fields.st.at("no3")    ->fld.data(), fields.sp.at("no3")->fld.data(), 
	    fields.st.at("no2")    ->fld.data(), fields.sp.at("no2")->fld.data(), 
	    fields.st.at("oh")     ->fld.data(), fields.sp.at("oh")->fld.data(), 
	    jval,emval, vdo3.data(),vdno.data(),vdno2.data(),vdhno3.data(),vdh2o2.data(),vdrooh.data(),vdhcho.data(),
	    rfa.data(), trfa,
	    fields.sp.at("qt") ->fld.data(),
	    Temp ->fld.data(), dt, sdt, switch_dt,
	    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
	    gd.icells, gd.ijcells, gd.dz.data(), fields.rhoref.data());
    fields.release_tmp(Temp);

//    isop_stat<TF>(
//	    fields.st.at("isop")   ->fld.data(), fields.sp.at("isop")->fld.data(), 
//	    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
//	    gd.icells, gd.ijcells);
}
#endif

template class Chemistry<double>;
//:template class Chemistry<float>;
