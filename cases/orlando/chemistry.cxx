
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


namespace 
{
#include "../cases/orlando/include/orlando_ohss_Parameters.h"
#include "../cases/orlando/include/orlando_ohss_Global.h"
#include "../cases/orlando/include/orlando_ohss_Sparse.h"
#include "../cases/orlando/include/orlando_ohss_Integrator.c"      /* needs to be modified */
#include "../cases/orlando/include/orlando_ohss_Function.c"        /* needs to be modified */
#include "../cases/orlando/include/orlando_ohss_LinearAlgebra.c"
#include "../cases/orlando/include/orlando_ohss_JacobianSP.c"
#include "../cases/orlando/include/orlando_ohss_Jacobian.c"
	

double C[NSPEC];                         /* Concentration of all species */
double * VAR = & C[0];
double * FIX = & C[18];
double RCONST[NREACT];                   /* Rate constants (global) */
double RF[NREACT];                       /* Modified: These contain the reaction fluxes also modify call in Integrator */
double Vdot[NVAR];                       /* Time derivative of variable species concentrations */
double TIME;                             /* Current integration time */
//double SUN;                              /* Sunlight intensity between [0,1] */
//double TEMP;                             /* Temperature */
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
	TF  ARRM( const TF A0, const TF B0, const TF C0, const TF TEMP )
	{
	return A0 * exp( -B0/TEMP ) 
		* pow( (TEMP/300.0), C0 );   
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
	TF usr_HO2_HO2(const  TF TEMP, const TF C_M, const TF C_H2O )
	/* for cesm-consistent reaction labels, equivalent to usr9 */
	/* HO2+HO2 -> H2O2+O2 */
	/* H2O included in fc calculation, not as reactant */
	{
	TF ko;
	TF kinf;
	TF fc;
	TF kr;

	if( C_H2O > (TF)0.0 ) 
	  {
	  ko   = (TF)2.3e-13 * exp( (TF)600./TEMP );
	  kinf = (TF)1.7e-33 * C_M * exp( (TF)1000./TEMP );
	  fc = (TF)1.0 + (TF)1.4e-21 *C_H2O* exp( (TF)2200./TEMP );
	  kr = (ko + kinf) * fc; 
	  }
	else
	  {
	  kr = (TF)0.0 ;
	  }
	return kr;
	}           

	template<typename TF>
	TF JPL_TROE( const TF k0_300K, const TF n, const TF kinf_300K, const TF m, 
		 const TF base, const TF temp, const TF cair )

	{
	/* !------------------------------------------------------------ */
	/* ! ... local variables */
	/* !------------------------------------------------------------ */
	TF zt_help;
	TF k0_T;
	TF kinf_T;
	TF k_ratio;

	zt_help = (TF)300./temp;
	k0_T    = k0_300K   * pow(zt_help,n) * cair ; /* ! k_0   at current T */
	kinf_T  = kinf_300K * pow(zt_help,m)        ; /* ! k_inf at current T */
	k_ratio = k0_T/kinf_T;

	return  k0_T/((TF)1.+k_ratio) * 
	      pow(base,(TF)1.0/
	      ((TF)1.0+pow(log10(k_ratio),(TF)2.0)));
	}           

	template<typename TF>
	TF usr_CO_OH_a( const TF temp, const TF c_m )
	/* ! for cesm-consistent reaction labels, equivalent to usr8 */
	/* ! CO+OH -> CO2+HO2 */
	{
	TF boltz = 1.38044e-16;
	return (TF)1.5e-13 * ((TF)1.+ (TF)6.e-7*boltz*c_m*temp);
	}

	template<typename TF>
	TF usr_N2O5_H2O( const TF k, const TF c_h2o )
	{
	return k * c_h2o;
	}

	template<typename TF>
	TF ARR2M( TF A0, TF B0, TF TEMP)
	{
	return A0 * exp( -B0/TEMP );
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
	    printf("%i  %12.3e ppb  %12.3e ppb/hour \n",k,iso_mean/cnt,iso_tend*3600.0/cnt);
	}
    }

    template<typename TF>
    void pss(
            TF* restrict th2o2, const TF* const restrict h2o2, 
            TF* restrict tch4, const TF* const restrict ch4, 
            TF* restrict tn2o5, const TF* const restrict n2o5, 
            TF* restrict thald, const TF* const restrict hald, 
            TF* restrict tco, const TF* const restrict co, 
            TF* restrict thcho, const TF* const restrict hcho, 
            TF* restrict tisopooh, const TF* const restrict isopooh, 
            TF* restrict tisop, const TF* const restrict isop, 
            TF* restrict tmvkmacr, const TF* const restrict mvkmacr, 
            TF* restrict txo2, const TF* const restrict xo2, 
            TF* restrict tisopao2, const TF* const restrict isopao2, 
            TF* restrict tno2, const TF* const restrict no2, 
            TF* restrict to3, const TF* const restrict o3, 
            TF* restrict tno, const TF* const restrict no, 
            TF* restrict tch3o2, const TF* const restrict ch3o2, 
            TF* restrict tisopbo2, const TF* const restrict isopbo2, 
            TF* restrict tno3, const TF* const restrict no3, 
            TF* restrict tho2, const TF* const restrict ho2, 
            TF* restrict oh,
	    const TF* const restrict jval, const TF* const restrict emval,
	    TF* restrict rfa, TF& trfa,
	    const TF* restrict qt,
	    const TF* restrict Temp, const TF rkdt, const TF switch_dt,
	    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk, const TF* const restrict dz, TF* const restrict rhoref)
    {
	    
	const int Pj_o31d = 0;
	const int Pj_h2o2 = 1;
	const int Pj_no2  = 2;
	const int Pj_no3a = 3;
	const int Pj_no3b = 4;
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
	TF vdo3   = (0.0);
	TF vdh2o2 = (0.0);
	TF vdno   = (0.0);
	TF vdno2  = (0.0);
	TF vdhcho = (0.0);
	TF vdisopooh = (0.0);
	TF vdhald    = (0.0);
	TF emvkmacr = (0.0);
	TF eno      = (0.0);

	for( int i = 0; i < NVAR; i++ ) {
	  RTOL[i] = RTOLS;
	  ATOL[i] = 1.0;
	  // tscale[i] = 1e20;
	}
	TF STEPMIN = 0.01;
	TF STEPMAX = 90;

	VAR = &C[0];
	FIX = &C[18];
        int nkpp = 0;
	int nderiv = 0;
	TF eisop;
	// update the time integration of the reaction fluxes:
	trfa += rkdt;
        for (int k=kstart; k<kend; ++k)
            {	
	    C_M = (TF)1e-3*rhoref[k]*Na/xmair;   // molecules/cm3 for chemistry!
	    const TF CFACTOR = C_M*(TF)1e-9 ;               // from ppb (units mixing ratio) to molecules/cm3
            if (k==kstart) {
	        vdo3   = TF(0.0056)/dz[k];   // 1/s
	        vdh2o2 = TF(0.0059)/dz[k];   // 1/s
	        vdno   = TF(0.0001)/dz[k];   // 1/s
	        vdno2  = TF(0.0027)/dz[k];   // 1/s
		vdhcho = TF(0.0032)/dz[k];   // 1/s
		vdisopooh = TF(0.0250)/dz[k];   // 1/s
		vdhald    = TF(0.0032)/dz[k];   // 1/s
		// emission/deposition fluxes:
		emvkmacr = TF(-0.004)*CFACTOR/dz[k]; // #/cm3 per s, seen orlando protocol
		eisop    = emval[0]*CFACTOR/dz[k];
		eno      = TF(0.03)*CFACTOR/dz[k];
	    }
            else {
	        vdo3   = TF(0.0);
	        vdh2o2 = TF(0.0);
	        vdno   = TF(0.0);
	        vdno2  = TF(0.0);
		vdhcho = TF(0.0);
		vdisopooh = TF(0.0);
		vdhald    = TF(0.0);
		emvkmacr= TF(0.0);
		eisop   = TF(0.0);
		eno     = TF(0.0);
	    }
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
		    const TF C_H2O = std::max(qt[ijk]*xmair*C_M/xmh2o,(TF)1.0);                   // kg/kg --> molH2O/molAir --*C_M--> molecules/cm3 limit to 1 molecule/cm3 to avoid error usr_HO2_HO2
		    const TF SUN = 1.0; 
		    const TF TEMP = Temp[ijk];
		    //const TF TEMP = 298.0;
		    // convert to molecules per cm3:
                    VAR[ind_H2O2]    = std::max(h2o2[ijk]*CFACTOR,(TF)0.0);
		    VAR[ind_CH4 ]    = std::max(ch4[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_N2O5]    = std::max(n2o5[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_HALD]    = std::max(hald[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_CO  ]    = std::max(co[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_HCHO]    = std::max(hcho[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_ISOPOOH] = std::max(isopooh[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_ISOP]    = std::max(isop[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_MVKMACR] = std::max(mvkmacr[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_XO2]     = std::max(xo2[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_ISOPAO2] = std::max(isopao2[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_NO2]     = std::max(no2[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_O3]      = std::max(o3[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_NO]      = std::max(no[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_CH3O2]   = std::max(ch3o2[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_ISOPBO2] = std::max(isopbo2[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_NO3]     = std::max(no3[ijk]*CFACTOR,(TF)0.0);
                    VAR[ind_HO2]     = std::max(ho2[ijk]*CFACTOR,(TF)0.0);
		    
                    RCONST[0] = (usr_O3_hv_H2O(TEMP,C_M,C_H2O,SUN*jval[Pj_o31d]));
                    RCONST[1] = (SUN*jval[Pj_no2]);
                    RCONST[2] = (SUN*jval[Pj_ch2or]);
                    RCONST[3] = (SUN*jval[Pj_ch2om]);
                    RCONST[4] = (ARR2M((TF)4.8e-11,-(TF)250.0,TEMP));
                    RCONST[5] = (ARR2M((TF)1.70e-12,(TF)940.0,TEMP));
                    RCONST[6] = (usr_HO2_HO2(TEMP,C_M,C_H2O));
                    RCONST[7] = (ARRM((TF)2.03e-16,-(TF)693.0,(TF)4.57,TEMP));
                    RCONST[8] = (ARR2M((TF)1.4e-12,(TF)1310.0,TEMP));
                    RCONST[9] = (JPL_TROE((TF)1.8e-30,(TF)3.,(TF)2.8e-11,(TF)0.0,(TF)0.6,TEMP,C_M));
                    RCONST[10] = (ARR2M((TF)1.4e-13,(TF)2470.0,TEMP));
                    RCONST[11] = (ARR2M((TF)3.50e-12,-(TF)270.0,TEMP));
                    RCONST[12] = (JPL_TROE((TF)2.4e-30,(TF)3.0,(TF)1.6e-12,-(TF)0.10,
                                (TF)0.60,TEMP,C_M));
                    RCONST[13] = (ARR2M((TF)1.8E-11,-(TF)110.0,TEMP));
                    RCONST[14] = ((TF)4.0e-12);
		    RCONST[15] = RCONST[12]*(TF)1.724e+26 * exp(-(TF)10840.0 / TEMP);
                    RCONST[16] = (usr_N2O5_H2O((TF)2.5e-22,C_H2O));
                    RCONST[17] = (SUN*jval[Pj_no3a]+SUN*jval[Pj_no3b]);
                    RCONST[18] = (ARR2M((TF)2.80e-12,-(TF)300.0,TEMP));
                    RCONST[19] = (ARR2M((TF)4.10e-13,-(TF)750.0,TEMP));
                    RCONST[20] = (usr_CO_OH_a(TEMP,C_M));
                    RCONST[21] = (ARR2M((TF)2.45e-12,(TF)1775.0,TEMP));
                    RCONST[22] = (ARR2M((TF)5.50e-12,-(TF)125.0,TEMP));
                    RCONST[23] = (ARR2M((TF)2.54e-11,-(TF)410.0,TEMP));
                    RCONST[24] = (ARR2M((TF)5.0e-13,-(TF)400.0,TEMP));
                    RCONST[25] = (ARR2M((TF)8.0e-13,-(TF)700.0,TEMP));
                    RCONST[26] = (ARR2M((TF)4.4e-12,-(TF)180.0,TEMP));
                    RCONST[27] = (ARR2M((TF)5.0e-13,-(TF)400.0,TEMP));
                    RCONST[28] = (ARR2M((TF)8.0e-13,-(TF)700.0,TEMP));
                    RCONST[29] = (ARR2M((TF)1.6e+09,(TF)8300.0,TEMP));
                    RCONST[30] = (ARR2M((TF)4.4e-12,-(TF)180.0,TEMP));
                    RCONST[31] = (ARR2M((TF)1.52e-11,-(TF)200.0,TEMP));
                    RCONST[32] = (SUN*jval[Pj_ch3o2h]);
                    RCONST[33] = (ARR2M((TF)1.86e-11,-(TF)175.0,TEMP));
                    RCONST[34] = ((TF)0.004*SUN*jval[Pj_no2]);
                    RCONST[35] = (ARR2M((TF)8.0e-13,-(TF)700.0,TEMP));
                    RCONST[36] = (ARR2M((TF)2.7e-12,-(TF)360.0,TEMP));
                    RCONST[37] = ((TF)2.4e-11);
                    RCONST[38] = (ARR2M((TF)1.03E-14,(TF)1995.0,TEMP));
                    RCONST[39] = (ARR2M((TF)3.15e-12,(TF)450.0,TEMP));
                    RCONST[40] = (SUN*jval[Pj_h2o2]);
		    RCONST[41] = (vdh2o2);
		    RCONST[42] = (vdo3);
		    RCONST[43] = (vdno);
		    RCONST[44] = (vdno2);
		    RCONST[45] = (vdhcho);
		    RCONST[46] = (vdisopooh);
		    RCONST[47] = (vdhald);
		    RCONST[48] = (eno);
		    RCONST[49] = (eisop);
		    RCONST[50] = (emvkmacr);

		    FIX[0] = ((TF)2.0*RCONST[0]*VAR[12] + RCONST[7]*VAR[12]*VAR[17] + RCONST[11]*VAR[13]*VAR[17] + RCONST[14]*VAR[16]*VAR[17] + (TF)2.0*RCONST[40]*VAR[0])/
		    	        ( RCONST[4]*VAR[17] +  RCONST[5]*VAR[12] +     RCONST[9]*VAR[11] +  RCONST[20]*VAR[4] + RCONST[21]*VAR[1] + 
		       	          RCONST[22]*VAR[5] +  RCONST[23]*VAR[7] + 0.4*RCONST[31]*VAR[6] +  RCONST[33]*VAR[3] + RCONST[37]*VAR[8]);
		    FIX[1] = C_H2O;
                    FIX[2] = C_M;
		    FIX[3] = (TF)1.0;    // species added to emit   
		    oh[ijk] = FIX[0]/CFACTOR;

		    Fun(VAR,FIX,RCONST,Vdot,RF);

		    // get statitics for reacion fluxes:
		    for (int l=0;l<NREACT;++l) rfa[(k-kstart)*NREACT+l] +=  RF[l]*rkdt;
		    // now to check reaction constants:
		    // for (int l=0;l<NREACT;++l) rfa[(k-kstart)*NREACT+l] +=  RCONST[l]*rkdt;

		    //for (int l=0; l<NVAR; ++l) printf (" %i %13.3e %13.3e %13.3e \n", l,VAR[l],Vdot[l],VAR[l]/ABS(Vdot[l]));
		    TF mint = (TF)1e20;
		    for (int l=0; l<NVAR; ++l)
			    if (ABS(Vdot[l]) > (TF)1e-5 && VAR[l]> (TF)1e-5) mint = std::min(mint,VAR[l]/ABS(Vdot[l])); 
		    //printf (" %i %13.3e \n ", ijk, mint);

		    if (mint < switch_dt)
		    {
			    nkpp += 1;

			    WCOPY(NVAR,VAR,1,VAR0,1);
			    INTEGRATE(  (TF)0.0 , rkdt );
			    
			    th2o2[ijk] +=    (VAR[ind_H2O2]-VAR0[ind_H2O2])/(rkdt*CFACTOR);
			    tch4[ijk] +=     (VAR[ind_CH4]-VAR0[ind_CH4])/(rkdt*CFACTOR);
			    tn2o5[ijk] +=    (VAR[ind_N2O5]-VAR0[ind_N2O5])/(rkdt*CFACTOR);
			    thald[ijk] +=    (VAR[ind_HALD]-VAR0[ind_HALD])/(rkdt*CFACTOR);
			    tco[ijk] +=      (VAR[ind_CO]-VAR0[ind_CO])/(rkdt*CFACTOR);
			    thcho[ijk] +=    (VAR[ind_HCHO]-VAR0[ind_HCHO])/(rkdt*CFACTOR);
			    tisopooh[ijk] += (VAR[ind_ISOPOOH]-VAR0[ind_ISOPOOH])/(rkdt*CFACTOR);
			    tisop[ijk] +=    (VAR[ind_ISOP]-VAR0[ind_ISOP])/(rkdt*CFACTOR);
			    tmvkmacr[ijk] += (VAR[ind_MVKMACR]-VAR0[ind_MVKMACR])/(rkdt*CFACTOR);
			    txo2[ijk] +=     (VAR[ind_XO2]-VAR0[ind_XO2])/(rkdt*CFACTOR);
			    tisopao2[ijk] += (VAR[ind_ISOPAO2]-VAR0[ind_ISOPAO2])/(rkdt*CFACTOR);
			    tno2[ijk] +=     (VAR[ind_NO2]-VAR0[ind_NO2])/(rkdt*CFACTOR);
			    to3[ijk] +=      (VAR[ind_O3]-VAR0[ind_O3])/(rkdt*CFACTOR);
			    tno[ijk] +=      (VAR[ind_NO]-VAR0[ind_NO])/(rkdt*CFACTOR);
			    tch3o2[ijk] +=   (VAR[ind_CH3O2]-VAR0[ind_CH3O2])/(rkdt*CFACTOR);
			    tisopbo2[ijk] += (VAR[ind_ISOPBO2]-VAR0[ind_ISOPBO2])/(rkdt*CFACTOR);
			    tno3[ijk] +=     (VAR[ind_NO3]-VAR0[ind_NO3])/(rkdt*CFACTOR);
			    ///  !!tho2[ijk] +=     (VAR[ind_HO2]-VAR0[ind_HO2])/(rkdt*CFACTOR);
			    tho2[ijk] +=     (VAR[ind_HO2]-VAR0[ind_HO2])/(rkdt*CFACTOR);
		    }
		    else
		    {
			    nderiv += 1;
			    th2o2[ijk] +=    Vdot[ind_H2O2]/CFACTOR;
   			    tch4[ijk] +=     Vdot[ind_CH4]/CFACTOR;
			    tn2o5[ijk] +=    Vdot[ind_N2O5]/CFACTOR;
			    thald[ijk] +=    Vdot[ind_HALD]/CFACTOR;
			    tco[ijk] +=      Vdot[ind_CO]/CFACTOR;
			    thcho[ijk] +=    Vdot[ind_HCHO]/CFACTOR;
			    tisopooh[ijk] += Vdot[ind_ISOPOOH]/CFACTOR;
			    tisop[ijk] +=    Vdot[ind_ISOP]/CFACTOR;
			    tmvkmacr[ijk] += Vdot[ind_MVKMACR]/CFACTOR;
			    txo2[ijk] +=     Vdot[ind_XO2]/CFACTOR;
			    tisopao2[ijk] += Vdot[ind_ISOPAO2]/CFACTOR;
			    tno2[ijk] +=     Vdot[ind_NO2]/CFACTOR;
			    to3[ijk] +=      Vdot[ind_O3]/CFACTOR;
			    tno[ijk] +=      Vdot[ind_NO]/CFACTOR;
			    tch3o2[ijk] +=   Vdot[ind_CH3O2]/CFACTOR;
			    tisopbo2[ijk] += Vdot[ind_ISOPBO2]/CFACTOR;
			    tno3[ijk] +=     Vdot[ind_NO3]/CFACTOR;
			    tho2[ijk] +=     Vdot[ind_HO2]/CFACTOR;
		    }
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
                }
	}
    }
}

template<typename TF>
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
	
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();
    fields.init_diagnostic_field("oh","oh","ppb", group_name, gd.sloc);
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

    stats.calc_stats("oh", *fields.sd.at("oh"), no_offset, no_threshold);

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
        else
            throw std::runtime_error("Invalid option for \"Chemistry type\"");
    }
    switch_dt = inputin.get_item<TF>("chemistry", "switch_dt","", (TF)1e5);
    statistics_counter = 0;

}

template <typename TF>
void Chemistry<TF>::create(const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    for (auto& it : fields.st) if( cmap[it.first].type == Chemistry_type::disabled) return;
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
    jno3a.resize(time_dim_length);
    jno3b.resize(time_dim_length);
    jch2or.resize(time_dim_length);
    jch2om.resize(time_dim_length);
    jch3o2h.resize(time_dim_length);
    emi_isop.resize(time_dim_length);

    group_nc.get_variable(time, time_dim, {0}, {time_dim_length});
    group_nc.get_variable(jo31d, jname[0],  {0}, {time_dim_length});
    group_nc.get_variable(jh2o2, jname[1],  {0}, {time_dim_length});
    group_nc.get_variable(jno2, jname[2],  {0}, {time_dim_length});
    group_nc.get_variable(jno3a, jname[3],  {0}, {time_dim_length});
    group_nc.get_variable(jno3b, jname[4],  {0}, {time_dim_length});
    group_nc.get_variable(jch2or, jname[5],  {0}, {time_dim_length});
    group_nc.get_variable(jch2om, jname[6],  {0}, {time_dim_length});
    group_nc.get_variable(jch3o2h, jname[7],  {0}, {time_dim_length});
    group_nc.get_variable(emi_isop, "emi_isop",  {0}, {time_dim_length});

    // Stats:
    const std::string group_name = "default";
    const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
    const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
    const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};

    // Add the profiles to te statistics
    if (stats.get_switch())
    {
        stats.add_profs(*fields.sd.at("oh"), "z", stat_op_w, group_name);
    }

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
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    for (auto& it : fields.st) if( cmap[it.first].type == Chemistry_type::disabled) return;

    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);
    jval[0] = ifac.fac0 * jo31d[ifac.index0] + ifac.fac1 * jo31d[ifac.index1];
    jval[1] = ifac.fac0 * jh2o2[ifac.index0] + ifac.fac1 * jh2o2[ifac.index1];
    jval[2] = ifac.fac0 * jno2[ifac.index0]  + ifac.fac1 * jno2[ifac.index1];
    jval[3] = ifac.fac0 * jno3a[ifac.index0] + ifac.fac1 * jno3a[ifac.index1];
    jval[4] = ifac.fac0 * jno3b[ifac.index0] + ifac.fac1 * jno3b[ifac.index1];
    jval[5] = ifac.fac0 * jch2or[ifac.index0] + ifac.fac1 * jch2or[ifac.index1];
    jval[6] = ifac.fac0 * jch2om[ifac.index0] + ifac.fac1 * jch2om[ifac.index1];
    jval[7] = ifac.fac0 * jch3o2h[ifac.index0] + ifac.fac1 * jch3o2h[ifac.index1];
    emval[0] = ifac.fac0 * emi_isop[ifac.index0] + ifac.fac1 * emi_isop[ifac.index1];
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

    // determine sub time step:
    TF rkdt = 0.0;
    if (abs(sdt/dt - 1./3.) < 1e-5)
	    rkdt = dt*(TF)1./(TF)3.;
    else if (abs(sdt/dt - 15./16.) < 1e-5)
	    rkdt = dt*(TF)5./(TF)12.;
    else if (abs(sdt/dt - 8./15.) < 1e-5)
	    rkdt = dt*(TF)1./(TF)4.;
    else
	    throw std::runtime_error("Invalid time step in RK3");
 
    pss<TF>(
	    fields.st.at("h2o2")   ->fld.data(), fields.sp.at("h2o2")->fld.data(), 
	    fields.st.at("ch4")    ->fld.data(), fields.sp.at("ch4")->fld.data(), 
	    fields.st.at("n2o5")   ->fld.data(), fields.sp.at("n2o5")->fld.data(), 
	    fields.st.at("hald")   ->fld.data(), fields.sp.at("hald")->fld.data(), 
	    fields.st.at("co")     ->fld.data(), fields.sp.at("co")->fld.data(), 
	    fields.st.at("hcho")   ->fld.data(), fields.sp.at("hcho")->fld.data(), 
	    fields.st.at("isopooh")->fld.data(), fields.sp.at("isopooh")->fld.data(), 
	    fields.st.at("isop")   ->fld.data(), fields.sp.at("isop")->fld.data(), 
	    fields.st.at("mvkmacr")->fld.data(), fields.sp.at("mvkmacr")->fld.data(), 
	    fields.st.at("xo2")    ->fld.data(), fields.sp.at("xo2")->fld.data(), 
	    fields.st.at("isopao2")->fld.data(), fields.sp.at("isopao2")->fld.data(), 
	    fields.st.at("no2")    ->fld.data(), fields.sp.at("no2")->fld.data(), 
	    fields.st.at("o3")     ->fld.data(), fields.sp.at("o3")->fld.data(), 
	    fields.st.at("no")     ->fld.data(), fields.sp.at("no")->fld.data(), 
	    fields.st.at("ch3o2")  ->fld.data(), fields.sp.at("ch3o2")->fld.data(), 
	    fields.st.at("isopbo2")->fld.data(), fields.sp.at("isopbo2")->fld.data(), 
	    fields.st.at("no3")    ->fld.data(), fields.sp.at("no3")->fld.data(), 
	    fields.st.at("ho2")    ->fld.data(), fields.sp.at("ho2")->fld.data(), 
	    fields.sd.at("oh")->fld.data(), 
	    jval,emval,
	    rfa.data(), trfa,
	    fields.st.at("qt") ->fld.data(),
	    Temp ->fld.data(), rkdt, switch_dt,
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
