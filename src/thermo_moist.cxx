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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <netcdf>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "diff_smag2.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"
#include "stats.h"
#include "master.h"
#include "cross.h"
#include "dump.h"
#include "thermo_moist_functions.h"
#include "timeloop.h"

using Finite_difference::O2::interp2;
using Finite_difference::O4::interp4;
using namespace Constants;
using namespace Thermo_moist_functions;

namespace
{
    // Settings for now here for convenience....
    const double pi      = std::acos(-1.);
    const double Nc0     = 70e6;     // Fixed cloud droplet number
    const double K_t     = 2.5e-2;   // Conductivity of heat [J/(sKm)]
    const double D_v     = 3.e-5;    // Diffusivity of water vapor [m2/s]
    const double rho_w   = 1.e3;     // Density water
    const double rho_0   = 1.225;    // SB06, p48
    const double pirhow  = pi * rho_w / 6.;
    const double mc_min  = 4.2e-15;  // Min mean mass of cloud droplet
    const double mc_max  = 2.6e-10;  // Max mean mass of cloud droplet
    const double mr_min  = mc_max;   // Min mean mass of precipitation drop
    const double mr_max  = 3e-6;     // Max mean mass of precipitation drop // as in UCLA-LES
    const double ql_min  = 1.e-6;    // Min cloud liquid water for which calculations are performed 
    const double qr_min  = 1.e-15;   // Min rain liquid water for which calculations are performed 

    // Rational tanh approximation  
    inline double tanh2(const double x)
    {
        return x * (27 + x * x) / (27 + 9 * x * x);
    }

    // Given rain water content (qr), number density (nr) and density (rho)
    // calculate mean mass of rain drop
    inline double calc_rain_mass(const double qr, const double nr, const double rho)
    {
        //double mr = rho * qr / (nr + dsmall);
        double mr = rho * qr / std::max(nr, 1.);
        mr        = std::min(std::max(mr, mr_min), mr_max);
        return mr;
    }

    // Given mean mass rain drop, calculate mean diameter
    inline double calc_rain_diameter(const double mr)
    {
        return pow(mr/pirhow, 1./3.);
    }

    // Shape parameter mu_r
    inline double calc_mu_r(const double dr)
    {
        //return 1./3.; // SB06
        //return 10. * (1. + tanh(1200 * (dr - 0.0015))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
        return 10. * (1. + tanh2(1200 * (dr - 0.0015))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
    }

    // Slope parameter lambda_r
    inline double calc_lambda_r(const double mur, const double dr)
    {
        return pow((mur+3)*(mur+2)*(mur+1), 1./3.) / dr;
    }

    inline double minmod(const double a, const double b)
    {
        return copysign(1., a) * std::max(0., std::min(std::abs(a), copysign(1., a)*b));
    }

    inline double min3(const double a, const double b, const double c)
    {
        return std::min(a,std::min(b,c));
    }

    inline double max3(const double a, const double b, const double c)
    {
        return std::max(a,std::max(b,c));
    }
}

// Microphysics per 2D (xz) slices
namespace mp2d
{
    // Calculate microphysics properties which are used in multiple routines
    void prepare_microphysics_slice(double* const restrict rain_mass, double* const restrict rain_diameter,
                                    double* const restrict mu_r, double* const restrict lambda_r,
                                    const double* const restrict qr, const double* const restrict nr,
                                    const double* const restrict rho, 
                                    const int istart, const int iend,
                                    const int kstart, const int kend,
                                    const int icells, const int ijcells, const int j)
    {
    
        for (int k=kstart; k<kend; k++) 
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                if(qr[ijk] > qr_min)
                {
                    rain_mass[ik]     = calc_rain_mass(qr[ijk], nr[ijk], rho[k]); 
                    rain_diameter[ik] = calc_rain_diameter(rain_mass[ik]);
                    mu_r[ik]          = calc_mu_r(rain_diameter[ik]);
                    lambda_r[ik]      = calc_lambda_r(mu_r[ik], rain_diameter[ik]);
                }
                else
                {
                    rain_mass[ik]     = 0; 
                    rain_diameter[ik] = 0;
                    mu_r[ik]          = 0;
                    lambda_r[ik]      = 0;
                }
            }
    }

    // Evaporation: evaporation of rain drops in unsaturated environment
    void evaporation(double* const restrict qrt, double* const restrict nrt,
                     double* const restrict qtt, double* const restrict thlt,
                     const double* const restrict qr, const double* const restrict nr,
                     const double* const restrict ql, const double* const restrict qt, const double* const restrict thl, 
                     const double* const restrict rho, const double* const restrict exner, const double* const restrict p,
                     const double* const restrict rain_mass, const double* const restrict rain_diameter,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int jj, const int kk, const int j)
    {
        const double lambda_evap = 1.; // 1.0 in UCLA, 0.7 in DALES

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*jj;
                const int ijk = i + j*jj + k*kk;

                if(qr[ijk] > qr_min)
                {
                    const double mr  = rain_mass[ik];
                    const double dr  = rain_diameter[ik];

                    const double T   = thl[ijk] * exner[k] + (Lv * ql[ijk]) / (cp * exner[k]); // Absolute temperature [K]
                    const double Glv = pow(Rv * T / (esat(T) * D_v) + (Lv / (K_t * T)) * (Lv / (Rv * T) - 1), -1); // Cond/evap rate (kg m-1 s-1)?

                    const double S   = (qt[ijk] - ql[ijk]) / qsat(p[k], T) - 1; // Saturation
                    const double F   = 1.; // Evaporation excludes ventilation term from SB06 (like UCLA, unimportant term? TODO: test)

                    const double ev_tend = 2. * pi * dr * Glv * S * F * nr[ijk] / rho[k];             

                    qrt[ijk]  += ev_tend;
                    nrt[ijk]  += lambda_evap * ev_tend * rho[k] / mr;
                    qtt[ijk]  -= ev_tend;
                    thlt[ijk] += Lv / (cp * exner[k]) * ev_tend; 
                }
            }
    }

    // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
    void selfcollection_breakup(double* const restrict nrt, const double* const restrict qr, const double* const restrict nr, const double* const restrict rho,
                                const double* const restrict rain_mass, const double* const restrict rain_diameter,
                                const double* const restrict lambda_r,
                                const int istart, const int jstart, const int kstart,
                                const int iend,   const int jend,   const int kend,
                                const int jj, const int kk, const int j)
    {
        const double k_rr     = 7.12;   // SB06, p49
        const double kappa_rr = 60.7;   // SB06, p49 
        const double D_eq     = 0.9e-3; // SB06, list of symbols 
        const double k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const double k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*jj;
                const int ijk = i + j*jj + k*kk;

                if(qr[ijk] > qr_min)
                {
                    // Calculate mean rain drop mass and diameter
                    const double mr      = rain_mass[ik];
                    const double dr      = rain_diameter[ik];
                    const double lambdar = lambda_r[ik];

                    // Selfcollection
                    const double sc_tend = -k_rr * nr[ijk] * qr[ijk]*rho[k] * pow(1. + kappa_rr / lambdar * pow(pirhow, 1./3.), -9) * pow(rho_0 / rho[k], 0.5);
                    nrt[ijk] += sc_tend; 

                    // Breakup
                    const double dDr = dr - D_eq;
                    if(dr > 0.35e-3)
                    {
                        double phi_br;
                        if(dr <= D_eq)
                            phi_br = k_br1 * dDr;
                        else
                            phi_br = 2. * exp(k_br2 * dDr) - 1.; 
                        const double br_tend = -(phi_br + 1.) * sc_tend;
                        nrt[ijk] += br_tend; 
                    }
                }
            }
    }

    // Sedimentation from Stevens and Seifert (2008)
    void sedimentation_ss08(double* const restrict qrt, double* const restrict nrt, 
                            double* const restrict tmpxz1, double* const restrict tmpxz2,
                            double* const restrict tmpxz3, double* const restrict tmpxz4,
                            double* const restrict tmpxz5, double* const restrict tmpxz6,
                            const double* const restrict mu_r, const double* const restrict lambda_r,
                            const double* const restrict qr, const double* const restrict nr, 
                            const double* const restrict rho, const double* const restrict dzi,
                            const double* const restrict dz, const double dt,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend,
                            const int icells, const int kcells, const int ijcells, const int j)
    {
        const double w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const double a_R = 9.65;   // SB06, p51
        const double c_R = 600;    // SB06, p51
        const double Dv  = 25.0e-6;
        const double b_R = a_R * exp(c_R*Dv); // UCLA-LES

        const int ikcells = icells * kcells;

        const int kk3d = ijcells;
        const int kk2d = icells;

        // 1. Calculate sedimentation velocity at cell center
        double* restrict w_qr = tmpxz1;
        double* restrict w_nr = tmpxz2;

        for (int k=kstart; k<kend; k++)
        {
            const double rho_n = pow(1.2 / rho[k], 0.5);
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;
              
                if(qr[ijk] > qr_min)
                {
                    // SS08:
                    w_qr[ik] = std::min(w_max, std::max(0.1, rho_n * a_R - b_R * pow(1. + c_R/lambda_r[ik], -1.*(mu_r[ik]+4))));
                    w_nr[ik] = std::min(w_max, std::max(0.1, rho_n * a_R - b_R * pow(1. + c_R/lambda_r[ik], -1.*(mu_r[ik]+1))));
                }
                else
                {
                    w_qr[ik] = 0.;
                    w_nr[ik] = 0.;
                }
            }
        }

        // 1.1 Set one ghost cell to zero
        for (int i=istart; i<iend; i++)
        {
            const int ik1  = i + (kstart-1)*icells;
            const int ik2  = i + (kend    )*icells;
            w_qr[ik1] = w_qr[ik1+kk2d];
            w_nr[ik1] = w_nr[ik1+kk2d];
            w_qr[ik2] = 0;
            w_nr[ik2] = 0;
        } 

         // 2. Calculate CFL number using interpolated sedimentation velocity
        double* restrict c_qr = tmpxz3;
        double* restrict c_nr = tmpxz4;

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*icells;
                c_qr[ik] = 0.25 * (w_qr[ik-kk2d] + 2.*w_qr[ik] + w_qr[ik+kk2d]) * dzi[k] * dt; 
                c_nr[ik] = 0.25 * (w_nr[ik-kk2d] + 2.*w_nr[ik] + w_nr[ik+kk2d]) * dzi[k] * dt; 
            }

        // 3. Calculate slopes
        double* restrict slope_qr = tmpxz1;
        double* restrict slope_nr = tmpxz2;

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                slope_qr[ik] = minmod(qr[ijk]-qr[ijk-kk3d], qr[ijk+kk3d]-qr[ijk]);
                slope_nr[ik] = minmod(nr[ijk]-nr[ijk-kk3d], nr[ijk+kk3d]-nr[ijk]);
            }

        // Calculate flux
        double* restrict flux_qr = tmpxz5;
        double* restrict flux_nr = tmpxz6;

        // Set the fluxes at the top of the domain (kend) to zero
        for (int i=istart; i<iend; i++)
        {
            const int ik  = i + kend*icells;
            flux_qr[ik] = 0;
            flux_nr[ik] = 0;
        } 

        for (int k=kend-1; k>kstart-1; k--)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {    
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                int kk;
                double ftot, dzz, cc;

                // q_rain
                kk    = k;  // current grid level
                ftot  = 0;  // cumulative 'flux' (kg m-2) 
                dzz   = 0;  // distance from zh[k]
                cc    = std::min(1., c_qr[ik]);
                while (cc > 0 && kk < kend)
                {
                    const int ikk  = i + kk*icells;
                    const int ijkk = i + j*icells + kk*ijcells;

                    ftot  += rho[kk] * (qr[ijkk] + 0.5 * slope_qr[ikk] * (1.-cc)) * cc * dz[kk];

                    dzz   += dz[kk];
                    kk    += 1;
                    cc     = std::min(1., c_qr[ikk] - dzz*dzi[kk]);
                }

                // Given flux at top, limit bottom flux such that the total rain content stays >= 0.
                ftot = std::min(ftot, rho[k] * dz[k] * qr[ijk] - flux_qr[ik+icells] * dt);
                flux_qr[ik] = -ftot / dt;

                // number density
                kk    = k;  // current grid level
                ftot  = 0;  // cumulative 'flux'
                dzz   = 0;  // distance from zh[k]
                cc    = std::min(1., c_nr[ik]);
                while (cc > 0 && kk < kend)
                {    
                    const int ikk  = i + kk*icells;
                    const int ijkk = i + j*icells + kk*ijcells;

                    ftot += rho[kk] * (nr[ijkk] + 0.5 * slope_nr[ikk] * (1.-cc)) * cc * dz[kk];

                    dzz   += dz[kk];
                    kk    += 1;
                    cc     = std::min(1., c_nr[ikk] - dzz*dzi[k]);
                }

                // Given flux at top, limit bottom flux such that the number density stays >= 0.
                ftot = std::min(ftot, rho[k] * dz[k] * nr[ijk] - flux_nr[ik+icells] * dt);
                flux_nr[ik] = -ftot / dt;
            }

        // Calculate tendency
        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {    
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                qrt[ijk] += -(flux_qr[ik+kk2d] - flux_qr[ik]) / rho[k] * dzi[k]; 
                nrt[ijk] += -(flux_nr[ik+kk2d] - flux_nr[ik]) / rho[k] * dzi[k]; 
            }
    }
}

namespace mp
{
    // Remove negative values from 3D field
    void remove_neg_values(double* const restrict field,
                           const int istart, const int jstart, const int kstart,
                           const int iend,   const int jend,   const int kend,
                           const int jj, const int kk)
    {
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    field[ijk] = std::max(0., field[ijk]);
                }
    }

    void zero(double* const restrict field, const int ncells)
    {
        for (int n=0; n<ncells; ++n)
            field[n] = 0.;
    }


    // Autoconversion: formation of rain drop by coagulating cloud droplets
    void autoconversion(double* const restrict qrt, double* const restrict nrt,
                        double* const restrict qtt, double* const restrict thlt,
                        const double* const restrict qr,  const double* const restrict ql, 
                        const double* const restrict rho, const double* const restrict exner,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int jj, const int kk)
    {
        const double x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES
        const double k_cc   = 9.44e9;        // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48
        const double nu_c   = 1;             // SB06, Table 1., same as UCLA-LES
        const double kccxs  = k_cc / (20. * x_star) * (nu_c+2)*(nu_c+4) / pow(nu_c+1, 2); 

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(ql[ijk] > ql_min)
                    {
                        const double xc      = rho[k] * ql[ijk] / Nc0; // Mean mass of cloud drops [kg]
                        const double tau     = 1. - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall); // SB06, Eq 5
                        const double phi_au  = 600. * pow(tau, 0.68) * pow(1. - pow(tau, 0.68), 3); // UCLA-LES
                        //const double phi_au  = 400. * pow(tau, 0.7) * pow(1. - pow(tau, 0.7), 3); // SB06, Eq 6
                        const double au_tend = rho[k] * kccxs * pow(ql[ijk], 2) * pow(xc, 2) *
                                               (1. + phi_au / pow(1-tau, 2)); // SB06, eq 4

                        qrt[ijk]  += au_tend; 
                        nrt[ijk]  += au_tend * rho[k] / x_star;
                        qtt[ijk]  -= au_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * au_tend; 
                    }
                }
    }

    // Evaporation: evaporation of rain drops in unsaturated environment
    void evaporation(double* const restrict qrt, double* const restrict nrt,
                     double* const restrict qtt, double* const restrict thlt,
                     const double* const restrict qr, const double* const restrict nr,
                     const double* const restrict ql, const double* const restrict qt, const double* const restrict thl, 
                     const double* const restrict rho, const double* const restrict exner, const double* const restrict p,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int jj, const int kk)
    {
        const double lambda_evap = 1.; // 1.0 in UCLA, 0.7 in DALES

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(qr[ijk] > qr_min)
                    {
                        // Calculate mean rain drop mass and diameter
                        const double mr  = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const double dr  = calc_rain_diameter(mr);

                        const double T   = thl[ijk] * exner[k] + (Lv * ql[ijk]) / (cp * exner[k]); // Absolute temperature [K]
                        const double Glv = pow(Rv * T / (esat(T) * D_v) + (Lv / (K_t * T)) * (Lv / (Rv * T) - 1), -1); // Cond/evap rate (kg m-1 s-1)?

                        const double S   = (qt[ijk] - ql[ijk]) / qsat(p[k], T) - 1; // Saturation
                        const double F   = 1.; // Evaporation excludes ventilation term from SB06 (like UCLA, unimportant term? TODO: test)

                        const double ev_tend = 2. * pi * dr * Glv * S * F * nr[ijk] / rho[k];             

                        qrt[ijk]  += ev_tend;
                        nrt[ijk]  += lambda_evap * ev_tend * rho[k] / mr;
                        qtt[ijk]  -= ev_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * ev_tend; 
                    }
                }
    }

    // Accreation: growth of raindrops collecting cloud droplets
    void accretion(double* const restrict qrt, double* const restrict qtt, double* const restrict thlt,
                   const double* const restrict qr,  const double* const restrict ql, 
                   const double* const restrict rho, const double* const restrict exner,
                   const int istart, const int jstart, const int kstart,
                   const int iend,   const int jend,   const int kend,
                   const int jj, const int kk)
    {
        const double k_cr  = 5.25; // SB06, p49

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(ql[ijk] > ql_min && qr[ijk] > qr_min)
                    {
                        const double tau     = 1 - ql[ijk] / (ql[ijk] + qr[ijk]); // SB06, Eq 5
                        const double phi_ac  = pow(tau / (tau + 5e-5), 4); // SB06, Eq 8
                        const double ac_tend = k_cr * ql[ijk] *  qr[ijk] * phi_ac * pow(rho_0 / rho[k], 0.5); // SB06, Eq 7 

                        qrt[ijk]  += ac_tend;
                        qtt[ijk]  -= ac_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * ac_tend; 
                    }
                }
    }


    // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
    void selfcollection_breakup(double* const restrict nrt, const double* const restrict qr, const double* const restrict nr, const double* const restrict rho,
                                const int istart, const int jstart, const int kstart,
                                const int iend,   const int jend,   const int kend,
                                const int jj, const int kk)
    {
        const double k_rr     = 7.12;   // SB06, p49
        const double kappa_rr = 60.7;   // SB06, p49 
        const double D_eq     = 0.9e-3; // SB06, list of symbols 
        const double k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const double k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(qr[ijk] > qr_min)
                    {
                        // Calculate mean rain drop mass and diameter
                        const double mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const double dr      = calc_rain_diameter(mr);
                        const double mur     = calc_mu_r(dr);
                        const double lambdar = calc_lambda_r(mur, dr);

                        // Selfcollection
                        const double sc_tend = -k_rr * nr[ijk] * qr[ijk]*rho[k] * pow(1. + kappa_rr / lambdar * pow(pirhow, 1./3.), -9) * pow(rho_0 / rho[k], 0.5);
                        nrt[ijk] += sc_tend; 

                        // Breakup
                        const double dDr = dr - D_eq;
                        if(dr > 0.35e-3)
                        {
                            double phi_br;
                            if(dr <= D_eq)
                                phi_br = k_br1 * dDr;
                            else
                                phi_br = 2. * exp(k_br2 * dDr) - 1.; 
                            const double br_tend = -(phi_br + 1.) * sc_tend;
                            nrt[ijk] += br_tend; 
                        }
                    }
                }
    }

    // Sedimentation from Stevens and Seifert (2008)
    void sedimentation_ss08(double* const restrict qrt, double* const restrict nrt, 
                            double* const restrict tmp1, double* const restrict tmp2,
                            const double* const restrict qr, const double* const restrict nr, 
                            const double* const restrict rho, const double* const restrict dzi,
                            const double* const restrict dz, const double dt,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend,
                            const int icells, const int kcells, const int ijcells)
    {
        const double w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const double a_R = 9.65;   // SB06, p51
        const double c_R = 600;    // SB06, p51
        const double Dv  = 25.0e-6;
        const double b_R = a_R * exp(c_R*Dv); // UCLA-LES

        const int ikcells = icells * kcells;

        for (int j=jstart; j<jend; ++j)
        {
            // 1. Calculate sedimentation velocity at cell center
            double* restrict w_qr = &tmp1[0*ikcells];
            double* restrict w_nr = &tmp1[1*ikcells];

            for (int k=kstart; k<kend; k++)
            {
                const double rho_n = pow(1.2 / rho[k], 0.5);
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;
                  
                    if(qr[ijk] > qr_min)
                    {
                        // Calculate mean rain drop mass and diameter
                        const double mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const double dr      = calc_rain_diameter(mr);
                        const double mur     = calc_mu_r(dr);
                        const double lambdar = calc_lambda_r(mur, dr);
              
                        // SS08:
                        w_qr[ik] = std::min(w_max, std::max(0.1, rho_n * a_R - b_R * pow(1. + c_R/lambdar, -1.*(mur+4))));
                        w_nr[ik] = std::min(w_max, std::max(0.1, rho_n * a_R - b_R * pow(1. + c_R/lambdar, -1.*(mur+1))));
                    }
                    else
                    {
                        w_qr[ik] = 0.;
                        w_nr[ik] = 0.;
                    }
                }
            }

            // 1.1 Set one ghost cell to zero
            for (int i=istart; i<iend; i++)
            {
                const int ik1  = i + (kstart-1)*icells;
                const int ik2  = i + (kend    )*icells;
                w_qr[ik1] = w_qr[ik1+icells];
                w_nr[ik1] = w_nr[ik1+icells];
                w_qr[ik2] = 0;
                w_nr[ik2] = 0;
            } 

             // 2. Calculate CFL number using interpolated sedimentation velocity
            double* restrict c_qr = &tmp1[2*ikcells];
            double* restrict c_nr = &tmp2[0*ikcells];

            for (int k=kstart; k<kend; k++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ik  = i + k*icells;
                    c_qr[ik] = 0.25 * (w_qr[ik-icells] + 2.*w_qr[ik] + w_qr[ik+icells]) * dzi[k] * dt; 
                    c_nr[ik] = 0.25 * (w_nr[ik-icells] + 2.*w_nr[ik] + w_nr[ik+icells]) * dzi[k] * dt; 
                }

            // 3. Calculate slopes with slope limiter to prevent forming new maxima/minima using slopes
            double* restrict slope_qr = &tmp1[0*ikcells];
            double* restrict slope_nr = &tmp1[1*ikcells];

            // Dissable slope limiter near surface, i.e. assume that e.g. qr[kstart]-qr[kstart-1] == qr[kstart+1]-qr[kstart] 
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*icells + kstart*ijcells;
                const int ik  = i + kstart*icells;

                slope_qr[ik] = qr[ijk+ijcells] - qr[ijk];
                slope_nr[ik] = nr[ijk+ijcells] - nr[ijk];
            }

            for (int k=kstart+1; k<kend; k++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;

                    slope_qr[ik] = minmod(qr[ijk]-qr[ijk-ijcells], qr[ijk+ijcells]-qr[ijk]);
                    slope_nr[ik] = minmod(nr[ijk]-nr[ijk-ijcells], nr[ijk+ijcells]-nr[ijk]);
                }

            // Calculate flux
            double* restrict flux_qr = &tmp2[1*ikcells];
            double* restrict flux_nr = &tmp2[2*ikcells];

            // Set the fluxes at the top of the domain (kend) to zero
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + kend*icells;
                flux_qr[ik] = 0;
                flux_nr[ik] = 0;
            } 

            for (int k=kend-1; k>kstart-1; k--)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {    
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;

                    int kk;
                    double ftot, dzz, cc;

                    // q_rain
                    kk    = k;  // current grid level
                    ftot  = 0;  // cumulative 'flux' (kg m-2) 
                    dzz   = 0;  // distance from zh[k]
                    cc    = std::min(1., c_qr[ik]);
                    while (cc > 0 && kk < kend)
                    {
                        const int ikk  = i + kk*icells;
                        const int ijkk = i + j*icells + kk*ijcells;
                        
                        ftot += rho[kk] * (qr[ijkk] + 0.5 * slope_qr[ikk] * (1.-cc)) * cc * dz[kk];
                        
                        dzz   += dz[kk];
                        kk    += 1;
                        cc     = std::min(1., c_qr[ikk] - dzz*dzi[kk]);
                    }

                    // Given flux at top, limit bottom flux such that the total rain content stays >= 0.
                    ftot = std::min(ftot, rho[k] * dz[k] * qr[ijk] - flux_qr[ik+icells] * dt - dsmall);
                    flux_qr[ik] = -ftot / dt;

                    // number density
                    kk    = k;  // current grid level
                    ftot  = 0;  // cumulative 'flux'
                    dzz   = 0;  // distance from zh[k]
                    cc    = std::min(1., c_nr[ik]);
                    while (cc > 0 && kk < kend)
                    {    
                        const int ikk  = i + kk*icells;
                        const int ijkk = i + j*icells + kk*ijcells;
            
                        ftot += rho[kk] * (nr[ijkk] + 0.5 * slope_nr[ikk] * (1.-cc)) * cc * dz[kk];

                        dzz   += dz[kk];
                        kk    += 1;
                        cc     = std::min(1., c_nr[ikk] - dzz*dzi[kk]);
                    }

                    // Given flux at top, limit bottom flux such that the number density stays >= 0.
                    ftot = std::min(ftot, rho[k] * dz[k] * nr[ijk] - flux_nr[ik+icells] * dt - dsmall);
                    flux_nr[ik] = -ftot / dt;
                }

            // Calculate tendency
            for (int k=kstart; k<kend; k++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {    
                    const int ijk = i + j*icells + k*ijcells;
                    const int ik  = i + k*icells;

                    qrt[ijk] += -(flux_qr[ik+icells] - flux_qr[ik]) / rho[k] * dzi[k]; 
                    nrt[ijk] += -(flux_nr[ik+icells] - flux_nr[ik]) / rho[k] * dzi[k]; 
                }
        }
    }

    // Calculate maximum sedimentation velocity
    double calc_max_sedimentation_cfl(double* const restrict w_qr,
                                      const double* const restrict qr, const double* const restrict nr, 
                                      const double* const restrict rho, const double* const restrict dzi,
                                      const double dt,
                                      const int istart, const int jstart, const int kstart,
                                      const int iend,   const int jend,   const int kend,
                                      const int icells, const int ijcells)
    {
        const double w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const double a_R = 9.65;   // SB06, p51
        const double c_R = 600;    // SB06, p51
        const double Dv  = 25.0e-6;
        const double b_R = a_R * exp(c_R*Dv); // UCLA-LES

        // Calculate sedimentation velocity at cell centre
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if(qr[ijk] > qr_min)
                    {
                        // Calculate mean rain drop mass and diameter
                        const double mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const double dr      = calc_rain_diameter(mr);
                        const double mur     = calc_mu_r(dr);
                        const double lambdar = calc_lambda_r(mur, dr);
                
                        w_qr[ijk] = std::min(w_max, std::max(0.1, a_R - b_R * pow(1. + c_R/lambdar, -1.*(mur+4))));
                    }
                    else
                    {
                        w_qr[ijk] = 0.;
                    }
                }

        // Calculate maximum CFL based on interpolated velocity
        double cfl_max = 1e-5;
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    
                    const double cfl_qr = 0.25 * (w_qr[ijk-ijcells] + 2.*w_qr[ijk] + w_qr[ijk+ijcells]) * dzi[k] * dt;
                    cfl_max = std::max(cfl_max, cfl_qr);
                }

        return cfl_max;
    }

} // End namespace


Thermo_moist::Thermo_moist(Model* modelin, Input* inputin) : Thermo(modelin, inputin)
{
    swthermo = "moist";

    thl0 = 0;
    qt0  = 0;

    thvref  = 0;
    thvrefh = 0;
    exnref  = 0;
    exnrefh = 0;
    pref    = 0;
    prefh   = 0;

    thvref_g  = 0;
    thvrefh_g = 0;
    exnref_g  = 0;
    exnrefh_g = 0;
    pref_g    = 0;
    prefh_g   = 0;

    int nerror = 0;

    // Option to overrule the prognostic variable
    nerror += inputin->get_item(&thvar, "thermo", "progvar", "", "thl");  // defaults to thl

    // BvS:micro Get microphysics switch, and init rain and number density
    nerror += inputin->get_item(&swmicro, "thermo", "swmicro", "", "0");
    if(swmicro == "2mom_warm")
    {
        #ifdef USECUDA
        master->print_error("swmicro = \"2mom_warm\" not (yet) implemented in CUDA\n");
        throw 1;
        #endif

        nerror += inputin->get_item(&swmicrobudget, "thermo", "swmicrobudget", "", "0");
        nerror += inputin->get_item(&swmicrobudget, "thermo", "swmicrobudget", "", "0");
        nerror += inputin->get_item(&cflmax_micro,  "thermo", "cflmax_micro",  "", 2.);

        // The microphysics requires three additional tmp fields
        const int n_tmp = 7;
        fields->set_minimum_tmp_fields(n_tmp);

        fields->init_prognostic_field("qr", "Rain water mixing ratio", "kg kg-1");
        fields->init_prognostic_field("nr", "Number density rain", "m-3");
    }    

    // Initialize the prognostic fields
    fields->init_prognostic_field(thvar, "Liquid water potential temperature", "K");
    fields->init_prognostic_field("qt", "Total water mixing ratio", "kg kg-1");

    // Get the diffusivities of temperature and moisture
    nerror += inputin->get_item(&fields->sp[thvar]->visc, "fields", "svisc", thvar );
    nerror += inputin->get_item(&fields->sp["qt"]->visc, "fields", "svisc", "qt");

    // Test if the diffusivities of theta and qt are equal, else throw error
    if (fields->sp[thvar]->visc != fields->sp["qt"]->visc)
    {
        master->print_error("The diffusivities of temperature and moisture must be equal\n");
        throw 1;
    }

    nerror += inputin->get_item(&pbot, "thermo", "pbot", "");

    // Get base state option (boussinesq or anelastic)
    nerror += inputin->get_item(&swbasestate, "thermo", "swbasestate", "", "");

    if (!(swbasestate == "boussinesq" || swbasestate == "anelastic"))
    {
        master->print_error("\"%s\" is an illegal value for swbasestate\n", swbasestate.c_str());
        throw 1;
    }

    // Note BvS: the base state calculation was never tested/supported for 4th order DNS,
    // and has therefor been removed for now.
    if (grid->swspatialorder == "4") // && swbasestate == "anelastic")
    {
        //master->print_error("Anelastic mode is not supported for swspatialorder=4\n");
        master->print_error("Moist thermodynamics are not supported for swspatialorder=4\n");
        throw 1;
    }

    // BvS test for updating hydrostatic prssure during run
    // swupdate..=0 -> initial base state pressure used in saturation calculation
    // swupdate..=1 -> base state pressure updated before saturation calculation
    nerror += inputin->get_item(&swupdatebasestate, "thermo", "swupdatebasestate", "");

    // Remove the data from the input that is not used, to avoid warnings.
    if (master->mode == "init")
    {
        inputin->flag_as_used("thermo", "thvref0");
        inputin->flag_as_used("thermo", "pbot");
    }

    if (nerror)
        throw 1;
}

Thermo_moist::~Thermo_moist()
{
    delete[] thl0;
    delete[] qt0;
    delete[] thvref;
    delete[] thvrefh;
    delete[] exnref;
    delete[] exnrefh;
    delete[] pref;
    delete[] prefh;

    #ifdef USECUDA
    clear_device();
    #endif
}

void Thermo_moist::init()
{
    stats = model->stats;

    thl0    = new double[grid->kcells];
    qt0     = new double[grid->kcells];
    thvref  = new double[grid->kcells];
    thvrefh = new double[grid->kcells];
    exnref  = new double[grid->kcells];
    exnrefh = new double[grid->kcells];
    pref    = new double[grid->kcells];
    prefh   = new double[grid->kcells];

    for (int k=0; k<grid->kcells; ++k)
    {
        thl0   [k] = 0.;
        qt0    [k] = 0.;
        thvref [k] = 0.;
        thvrefh[k] = 0.;
        exnref [k] = 0.;
        exnrefh[k] = 0.;
        pref   [k] = 0.;
        prefh  [k] = 0.;
    }

    init_cross();
    init_dump();
}

void Thermo_moist::create(Input* inputin)
{
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    // Enable automated calculation of horizontally averaged fields
    if (swupdatebasestate)
        fields->set_calc_mean_profs(true);

    // Calculate the base state profiles. With swupdatebasestate=1, these profiles are updated on every iteration.
    // 1. Take the initial profile as the reference
    if (inputin->get_prof(&thl0[grid->kstart], thvar, grid->kmax))
        throw 1;
    if (inputin->get_prof(&qt0[grid->kstart], "qt", grid->kmax))
        throw 1;

    // 2. Calculate surface and model top values thl and qt
    double thl0s, qt0s, thl0t, qt0t;
    thl0s = thl0[kstart] - grid->z[kstart]*(thl0[kstart+1]-thl0[kstart])*grid->dzhi[kstart+1];
    qt0s  = qt0[kstart]  - grid->z[kstart]*(qt0[kstart+1] -qt0[kstart] )*grid->dzhi[kstart+1];
    thl0t = thl0[kend-1] + (grid->zh[kend]-grid->z[kend-1])*(thl0[kend-1]-thl0[kend-2])*grid->dzhi[kend-1];
    qt0t  = qt0[kend-1]  + (grid->zh[kend]-grid->z[kend-1])*(qt0[kend-1]- qt0[kend-2] )*grid->dzhi[kend-1];

    // 3. Set the ghost cells for the reference temperature and moisture
    thl0[kstart-1]  = 2.*thl0s - thl0[kstart];
    thl0[kend]      = 2.*thl0t - thl0[kend-1];
    qt0[kstart-1]   = 2.*qt0s  - qt0[kstart];
    qt0[kend]       = 2.*qt0t  - qt0[kend-1];

    // 4. Calculate the initial/reference base state
    calc_base_state(pref, prefh, fields->rhoref, fields->rhorefh, thvref, thvrefh, exnref, exnrefh, thl0, qt0);

    // 5. In Boussinesq mode, overwrite reference temperature and density
    if (swbasestate == "boussinesq")
    {
        if (inputin->get_item(&thvref0, "thermo", "thvref0", ""))
            throw 1;

        for (int k=0; k<grid->kcells; ++k)
        {
            fields->rhoref[k]  = 1.;
            fields->rhorefh[k] = 1.;
            thvref[k]          = thvref0;
            thvrefh[k]         = thvref0;
        }
    }

    init_stat();
}

#ifndef USECUDA
void Thermo_moist::exec()
{
    const int kk = grid->ijcells;
    const int kcells = grid->kcells;

    // Re-calculate hydrostatic pressure and exner, pass dummy as rhoref,thvref to prevent overwriting base state
    double *tmp2 = fields->atmp["tmp2"]->data;
    if (swupdatebasestate)
        calc_base_state(pref, prefh,
                        &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells],
                        exnref, exnrefh, fields->sp[thvar]->datamean, fields->sp["qt"]->datamean);

    // extend later for gravity vector not normal to surface
    if (grid->swspatialorder == "2")
    {
        calc_buoyancy_tend_2nd(fields->wt->data, fields->sp[thvar]->data, fields->sp["qt"]->data, prefh,
                               &fields->atmp["tmp2"]->data[0*kk], &fields->atmp["tmp2"]->data[1*kk],
                               &fields->atmp["tmp2"]->data[2*kk], thvrefh);
    }
    //else if (grid->swspatialorder == "4")
    //{
    //    calc_buoyancy_tend_4th(fields->wt->data, fields->sp[thvar]->data, fields->sp["qt"]->data, prefh,
    //                           &fields->atmp["tmp2"]->data[0*kk], &fields->atmp["tmp2"]->data[1*kk],
    //                           &fields->atmp["tmp2"]->data[2*kk],
    //                           thvrefh);
    //}

    // 2-moment warm microphysics 
    if(swmicro == "2mom_warm")
        exec_microphysics();
}
#endif

unsigned long Thermo_moist::get_time_limit(unsigned long idt, const double dt)
{
    if(swmicro == "2mom_warm")
    {
        double cfl = mp::calc_max_sedimentation_cfl(fields->atmp["tmp1"]->data, fields->sp["qr"]->data, fields->sp["nr"]->data,
                                                    fields->rhoref, grid->dzi, dt,
                                                    grid->istart, grid->jstart, grid->kstart,
                                                    grid->iend,   grid->jend,   grid->kend,
                                                    grid->icells, grid->ijcells);
        grid->get_max(&cfl);
        return idt * cflmax_micro / cfl;
    }
    else
    {
        return Constants::ulhuge;
    }
}

// BvS:micro 
void Thermo_moist::exec_microphysics()
{
    // Switch to solve certain routines over xz-slices, to reduce calculations recurring in several microphysics routines
    bool per_slice = true;

    // Remove the negative values from the precipitation fields
    mp::remove_neg_values(fields->sp["qr"]->data, grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend, grid->icells, grid->ijcells);
    mp::remove_neg_values(fields->sp["nr"]->data, grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend, grid->icells, grid->ijcells);

    // Calculate the cloud liquid water concent using the saturation adjustment method
    calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp["thl"]->data, fields->sp["qt"]->data, pref);

    const double dt = model->timeloop->get_dt();

    // xz tmp slices for quantities which are used by multiple microphysics routines
    const int ikslice = grid->icells * grid->kcells;
    double* rain_mass = &fields->atmp["tmp2"]->data[0*ikslice]; 
    double* rain_diam = &fields->atmp["tmp2"]->data[1*ikslice]; 
    double* mu_r      = &fields->atmp["tmp2"]->data[2*ikslice]; 
    double* lambda_r  = &fields->atmp["tmp3"]->data[0*ikslice]; 

    // xz tmp slices for intermediate calculations
    double* tmpxz1    = &fields->atmp["tmp3"]->data[1*ikslice];
    double* tmpxz2    = &fields->atmp["tmp3"]->data[2*ikslice];
    double* tmpxz3    = &fields->atmp["tmp4"]->data[0*ikslice];
    double* tmpxz4    = &fields->atmp["tmp4"]->data[1*ikslice];
    double* tmpxz5    = &fields->atmp["tmp4"]->data[2*ikslice];
    double* tmpxz6    = &fields->atmp["tmp5"]->data[0*ikslice];

    // Autoconversion; formation of rain drop by coagulating cloud droplets
    mp::autoconversion(fields->st["qr"]->data, fields->st["nr"]->data, fields->st["qt"]->data, fields->st["thl"]->data,
                       fields->sp["qr"]->data, fields->atmp["tmp1"]->data, fields->rhoref, exnref,
                       grid->istart, grid->jstart, grid->kstart, 
                       grid->iend,   grid->jend,   grid->kend, 
                       grid->icells, grid->ijcells);

    // Accretion; growth of raindrops collecting cloud droplets
    mp::accretion(fields->st["qr"]->data, fields->st["qt"]->data, fields->st["thl"]->data,
                  fields->sp["qr"]->data, fields->atmp["tmp1"]->data, fields->rhoref, exnref,
                  grid->istart, grid->jstart, grid->kstart, 
                  grid->iend,   grid->jend,   grid->kend, 
                  grid->icells, grid->ijcells);

    if(per_slice)
    {
        for (int j=grid->jstart; j<grid->jend; ++j)
        {
            mp2d::prepare_microphysics_slice(rain_mass, rain_diam, mu_r, lambda_r, fields->sp["qr"]->data, fields->sp["nr"]->data, fields->rhoref,
                                             grid->istart, grid->iend, grid->kstart, grid->kend, grid->icells, grid->ijcells, j);

            // Evaporation; evaporation of rain drops in unsaturated environment
            mp2d::evaporation(fields->st["qr"]->data, fields->st["nr"]->data,  fields->st["qt"]->data, fields->st["thl"]->data,
                              fields->sp["qr"]->data, fields->sp["nr"]->data,  fields->atmp["tmp1"]->data,
                              fields->sp["qt"]->data, fields->sp["thl"]->data, fields->rhoref, exnref, pref,
                              rain_mass, rain_diam,
                              grid->istart, grid->jstart, grid->kstart, 
                              grid->iend,   grid->jend,   grid->kend, 
                              grid->icells, grid->ijcells, j);

            // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
            mp2d::selfcollection_breakup(fields->st["nr"]->data, fields->sp["qr"]->data, fields->sp["nr"]->data, fields->rhoref,
                                         rain_mass, rain_diam, lambda_r,
                                         grid->istart, grid->jstart, grid->kstart, 
                                         grid->iend,   grid->jend,   grid->kend, 
                                         grid->icells, grid->ijcells, j);

            // Sedimentation; sub-grid sedimentation of rain 
            mp2d::sedimentation_ss08(fields->st["qr"]->data, fields->st["nr"]->data, 
                                     tmpxz1, tmpxz2, tmpxz3, tmpxz4, tmpxz5, tmpxz6, mu_r, lambda_r,
                                     fields->sp["qr"]->data, fields->sp["nr"]->data, 
                                     fields->rhoref, grid->dzi, grid->dz, dt,
                                     grid->istart, grid->jstart, grid->kstart, 
                                     grid->iend,   grid->jend,   grid->kend, 
                                     grid->icells, grid->kcells, grid->ijcells, j);
        }
    }
    else
    {
        // Evaporation; evaporation of rain drops in unsaturated environment
        mp::evaporation(fields->st["qr"]->data, fields->st["nr"]->data,  fields->st["qt"]->data, fields->st["thl"]->data,
                        fields->sp["qr"]->data, fields->sp["nr"]->data,  fields->atmp["tmp1"]->data,
                        fields->sp["qt"]->data, fields->sp["thl"]->data, fields->rhoref, exnref, pref,
                        grid->istart, grid->jstart, grid->kstart, 
                        grid->iend,   grid->jend,   grid->kend, 
                        grid->icells, grid->ijcells);
       
        // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        mp::selfcollection_breakup(fields->st["nr"]->data, fields->sp["qr"]->data, fields->sp["nr"]->data, fields->rhoref,
                                   grid->istart, grid->jstart, grid->kstart, 
                                   grid->iend,   grid->jend,   grid->kend, 
                                   grid->icells, grid->ijcells);
    
        // Sedimentation; sub-grid sedimentation of rain 
        mp::sedimentation_ss08(fields->st["qr"]->data, fields->st["nr"]->data, 
                               fields->atmp["tmp4"]->data, fields->atmp["tmp5"]->data,
                               fields->sp["qr"]->data, fields->sp["nr"]->data, 
                               fields->rhoref, grid->dzi, grid->dz, dt,
                               grid->istart, grid->jstart, grid->kstart, 
                               grid->iend,   grid->jend,   grid->kend, 
                               grid->icells, grid->kcells, grid->ijcells);
    }
}

void Thermo_moist::get_mask(Field3d *mfield, Field3d *mfieldh, Mask *m)
{
    if (m->name == "ql")
    {
        calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
        calc_mask_ql(mfield->data, mfieldh->data, mfieldh->databot,
                     stats->nmask, stats->nmaskh, &stats->nmaskbot,
                     fields->atmp["tmp1"]->data);
    }
    else if (m->name == "qlcore")
    {
        calc_buoyancy(fields->atmp["tmp2"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp1"]->data,thvref);
        grid->calc_mean(fields->atmp["tmp2"]->datamean, fields->atmp["tmp2"]->data, grid->kcells);

        calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
        calc_mask_qlcore(mfield->data, mfieldh->data, mfieldh->databot,
                         stats->nmask, stats->nmaskh, &stats->nmaskbot,
                         fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp2"]->datamean);
    }
}

void Thermo_moist::calc_mask_ql(double* restrict mask, double* restrict maskh, double* restrict maskbot,
                                int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot,
                                double* restrict ql)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int k=grid->kstart; k<grid->kend; k++)
    {
        nmask[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = ql[ijk] > 0.;
                nmask[k] += ntmp;
                mask[ijk] = (double)ntmp;
            }
    }

    for (int k=grid->kstart; k<grid->kend+1; k++)
    {
        nmaskh[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = (ql[ijk-kk] + ql[ijk]) > 0.;

                nmaskh[k] += ntmp;
                maskh[ijk] = (double)ntmp;
            }
    }

    // Set the mask for surface projected quantities
    // In this case: ql at surface
    for (int j=grid->jstart; j<grid->jend; j++)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            maskbot[ij] = maskh[ijk];
        }

    grid->boundary_cyclic(mask);
    grid->boundary_cyclic(maskh);
    grid->boundary_cyclic_2d(maskbot);

    master->sum(nmask , grid->kcells);
    master->sum(nmaskh, grid->kcells);
    *nmaskbot = nmaskh[grid->kstart];

    // BvS: should no longer be necessary now that the ql ghost cells are set to zero
    //nmaskh[kstart] = 0;
    //nmaskh[kend  ] = 0;
}

void Thermo_moist::calc_mask_qlcore(double* restrict mask, double* restrict maskh, double* restrict maskbot,
                                    int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot,
                                    double* restrict ql, double* restrict b, double* restrict bmean)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int k=grid->kstart; k<grid->kend; k++)
    {
        nmask[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = (ql[ijk] > 0.)*(b[ijk]-bmean[k] > 0.);
                nmask[k] += ntmp;
                mask[ijk] = (double)ntmp;
            }
    }

    for (int k=grid->kstart; k<grid->kend+1; k++)
    {
        nmaskh[k] = 0;
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ntmp = (ql[ijk-kk]+ql[ijk] > 0.)*(b[ijk-kk]+b[ijk]-bmean[k-1]-bmean[k] > 0.);
                nmaskh[k] += ntmp;
                maskh[ijk] = (double)ntmp;
            }
    }

    // Set the mask for surface projected quantities
    // In this case: qlcore at surface
    for (int j=grid->jstart; j<grid->jend; j++)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            maskbot[ij] = maskh[ijk];
        }

    grid->boundary_cyclic(mask);
    grid->boundary_cyclic(maskh);
    grid->boundary_cyclic_2d(maskbot);

    master->sum(nmask , grid->kcells);
    master->sum(nmaskh, grid->kcells);
    *nmaskbot = nmaskh[grid->kstart];
}

void Thermo_moist::exec_stats(Mask *m)
{
    const double NoOffset = 0.;

    // calc the buoyancy and its surface flux for the profiles
    calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
    calc_buoyancy_fluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);

    // define location
    const int sloc[] = {0,0,0};

    // mean
    stats->calc_mean(m->profs["b"].data, fields->atmp["tmp1"]->data, NoOffset, sloc,
                     fields->atmp["tmp3"]->data, stats->nmask);

    // moments
    for (int n=2; n<5; ++n)
    {
        std::stringstream ss;
        ss << n;
        std::string sn = ss.str();
        stats->calc_moment(fields->atmp["tmp1"]->data, m->profs["b"].data, m->profs["b"+sn].data, n, sloc,
                           fields->atmp["tmp3"]->data, stats->nmask);
    }

    // calculate the gradients
    if (grid->swspatialorder == "2")
        stats->calc_grad_2nd(fields->atmp["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    else if (grid->swspatialorder == "4")
        stats->calc_grad_4th(fields->atmp["tmp1"]->data, m->profs["bgrad"].data, grid->dzhi4, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);

    // calculate turbulent fluxes
    if (grid->swspatialorder == "2")
        stats->calc_flux_2nd(fields->atmp["tmp1"]->data, m->profs["b"].data, fields->w->data, m->profs["w"].data,
                             m->profs["bw"].data, fields->atmp["tmp2"]->data, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    else if (grid->swspatialorder == "4")
        stats->calc_flux_4th(fields->atmp["tmp1"]->data, fields->w->data, m->profs["bw"].data, fields->atmp["tmp2"]->data, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);

    // calculate diffusive fluxes
    if (grid->swspatialorder == "2")
    {
        if (model->diff->get_switch() == "smag2")
        {
            Diff_smag_2 *diffptr = static_cast<Diff_smag_2 *>(model->diff);
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, fields->w->data, fields->sd["evisc"]->data,
                                 m->profs["bdiff"].data, grid->dzhi,
                                 fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->datafluxtop, diffptr->tPr, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
        }
        else
        {
            stats->calc_diff_2nd(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi, fields->sp[thvar]->visc, sloc,
                                 fields->atmp["tmp4"]->data, stats->nmaskh);
        }
    }
    else if (grid->swspatialorder == "4")
    {
        // take the diffusivity of temperature for that of buoyancy
        stats->calc_diff_4th(fields->atmp["tmp1"]->data, m->profs["bdiff"].data, grid->dzhi4, fields->sp[thvar]->visc, sloc,
                             fields->atmp["tmp4"]->data, stats->nmaskh);
    }

    // calculate the total fluxes
    stats->add_fluxes(m->profs["bflux"].data, m->profs["bw"].data, m->profs["bdiff"].data);

    // calculate the liquid water stats
    calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
    stats->calc_mean(m->profs["ql"].data, fields->atmp["tmp1"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
    stats->calc_count(fields->atmp["tmp1"]->data, m->profs["cfrac"].data, 0.,
                      fields->atmp["tmp3"]->data, stats->nmask);

    stats->calc_cover(fields->atmp["tmp1"]->data, fields->atmp["tmp4"]->databot, &stats->nmaskbot, &m->tseries["ccover"].data, 0.);
    stats->calc_path (fields->atmp["tmp1"]->data, fields->atmp["tmp4"]->databot, &stats->nmaskbot, &m->tseries["lwp"].data);

    // BvS:micro 
    if(swmicro == "2mom_warm")
    {
        stats->calc_path (fields->sp["qr"]->data, fields->atmp["tmp4"]->databot, &stats->nmaskbot, &m->tseries["rwp"].data);

        if(swmicrobudget == "1")
        {
            // Autoconversion
            mp::zero(fields->atmp["tmp2"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp5"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp6"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp7"]->data, grid->ncells);

            mp::autoconversion(fields->atmp["tmp2"]->data, fields->atmp["tmp5"]->data, fields->atmp["tmp6"]->data, fields->atmp["tmp7"]->data,
                               fields->sp["qr"]->data, fields->atmp["tmp1"]->data, fields->rhoref, exnref,
                               grid->istart, grid->jstart, grid->kstart, 
                               grid->iend,   grid->jend,   grid->kend, 
                               grid->icells, grid->ijcells);

            stats->calc_mean(m->profs["auto_qrt" ].data, fields->atmp["tmp2"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["auto_nrt" ].data, fields->atmp["tmp5"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["auto_qtt" ].data, fields->atmp["tmp6"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["auto_thlt"].data, fields->atmp["tmp7"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);

            // Evaporation
            mp::zero(fields->atmp["tmp2"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp5"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp6"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp7"]->data, grid->ncells);

            mp::evaporation(fields->atmp["tmp2"]->data, fields->atmp["tmp5"]->data,  fields->atmp["tmp6"]->data, fields->atmp["tmp7"]->data,
                            fields->sp["qr"]->data, fields->sp["nr"]->data,  fields->atmp["tmp1"]->data,
                            fields->sp["qt"]->data, fields->sp["thl"]->data, fields->rhoref, exnref, pref,
                            grid->istart, grid->jstart, grid->kstart, 
                            grid->iend,   grid->jend,   grid->kend, 
                            grid->icells, grid->ijcells);

            stats->calc_mean(m->profs["evap_qrt" ].data, fields->atmp["tmp2"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["evap_nrt" ].data, fields->atmp["tmp5"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["evap_qtt" ].data, fields->atmp["tmp6"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["evap_thlt"].data, fields->atmp["tmp7"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);

            // Accretion
            mp::zero(fields->atmp["tmp2"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp5"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp6"]->data, grid->ncells);

            mp::accretion(fields->atmp["tmp2"]->data, fields->atmp["tmp5"]->data, fields->atmp["tmp6"]->data,
                          fields->sp["qr"]->data, fields->atmp["tmp1"]->data, fields->rhoref, exnref,
                          grid->istart, grid->jstart, grid->kstart, 
                          grid->iend,   grid->jend,   grid->kend, 
                          grid->icells, grid->ijcells);

            stats->calc_mean(m->profs["accr_qrt" ].data, fields->atmp["tmp2"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["accr_qtt" ].data, fields->atmp["tmp5"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["accr_thlt"].data, fields->atmp["tmp6"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);

            // Selfcollection and breakup
            mp::zero(fields->atmp["tmp2"]->data, grid->ncells);

            mp::selfcollection_breakup(fields->atmp["tmp2"]->data, fields->sp["qr"]->data, fields->sp["nr"]->data, fields->rhoref,
                                       grid->istart, grid->jstart, grid->kstart, 
                                       grid->iend,   grid->jend,   grid->kend, 
                                       grid->icells, grid->ijcells);

            stats->calc_mean(m->profs["scbr_nrt" ].data, fields->atmp["tmp2"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);

            // Sedimentation
            mp::zero(fields->atmp["tmp2"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp5"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp6"]->data, grid->ncells);
            mp::zero(fields->atmp["tmp7"]->data, grid->ncells);

            // 1. Get number of substeps based on sedimentation with CFL=1
            //const double dt = model->timeloop->get_sub_time_step();
            //int nsubstep = mp::get_sedimentation_steps(fields->sp["qr"]->data, fields->sp["nr"]->data, 
            //                                           fields->rhoref, grid->dzi, dt,
            //                                           grid->istart, grid->jstart, grid->kstart, 
            //                                           grid->iend,   grid->jend,   grid->kend, 
            //                                           grid->icells, grid->ijcells);

            //// 2. Synchronize over all MPI processes 
            //grid->get_max(&nsubstep);

            //// Stay a bit further from CFL=1
            //nsubstep *= 2;

            //// Sedimentation in nsubstep steps:
            //mp::sedimentation_sub(fields->atmp["tmp2"]->data, fields->atmp["tmp5"]->data, 
            //                      fields->atmp["tmp6"]->data, fields->atmp["tmp7"]->data,
            //                      fields->sp["qr"]->data, fields->sp["nr"]->data, 
            //                      fields->rhoref, grid->dzi, grid->dzhi, dt,
            //                      grid->istart, grid->jstart, grid->kstart, 
            //                      grid->iend,   grid->jend,   grid->kend, 
            //                      grid->icells, grid->kcells, grid->ijcells,
            //                      nsubstep);

            const double dt = model->timeloop->get_sub_time_step();
            mp::sedimentation_ss08(fields->atmp["tmp2"]->data, fields->atmp["tmp5"]->data, 
                                   fields->atmp["tmp6"]->data, fields->atmp["tmp7"]->data,
                                   fields->sp["qr"]->data, fields->sp["nr"]->data, 
                                   fields->rhoref, grid->dzi, grid->dz, dt,
                                   grid->istart, grid->jstart, grid->kstart, 
                                   grid->iend,   grid->jend,   grid->kend, 
                                   grid->icells, grid->kcells, grid->ijcells);

            stats->calc_mean(m->profs["sed_qrt"].data, fields->atmp["tmp2"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
            stats->calc_mean(m->profs["sed_nrt"].data, fields->atmp["tmp5"]->data, NoOffset, sloc, fields->atmp["tmp3"]->data, stats->nmask);
        }
    }

    // Calculate base state in tmp array
    if (swupdatebasestate == 1)
    {
        const int kcells = grid->kcells;
        double* restrict tmp2 = fields->atmp["tmp2"]->data;
        calc_base_state(&tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells],
                        &tmp2[4*kcells], &tmp2[5*kcells], &tmp2[6*kcells], &tmp2[7*kcells],
                        fields->sp[thvar]->datamean, fields->sp["qt"]->datamean);

        for (int k=0; k<kcells; ++k)
        {
            m->profs["ph"  ].data[k] = tmp2[0*kcells+k];
            m->profs["phh" ].data[k] = tmp2[1*kcells+k];
            m->profs["rho" ].data[k] = tmp2[2*kcells+k];
            m->profs["rhoh"].data[k] = tmp2[3*kcells+k];
        }
    }
}

void Thermo_moist::exec_cross()
{
    int nerror = 0;

    Cross* cross = model->cross;

    // With one additional temp field, we wouldn't have to re-calculate the ql or b field for simple,lngrad,path, etc.
    for (std::vector<std::string>::iterator it=crosslist.begin(); it<crosslist.end(); ++it)
    {
        /* BvS: for now, don't call getThermoField() or getBuoyancySurf(), but directly the function itself. With CUDA enabled,
           statistics etc. is done on the host, while getThermoField() is executed on the GPU */

        if (*it == "b")
        {
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
            nerror += cross->cross_simple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
        }
        else if (*it == "ql")
        {
            calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
            nerror += cross->cross_simple(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, *it);
        }
        else if (*it == "blngrad")
        {
            calc_buoyancy(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp2"]->data, thvref);
            // Note: tmp1 twice used as argument -> overwritten in crosspath()
            nerror += cross->cross_lngrad(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, grid->dzi4, *it);
        }
        else if (*it == "qlpath")
        {
            calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
            // Note: tmp1 twice used as argument -> overwritten in crosspath()
            nerror += cross->cross_path(fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, "qlpath");
        }
        else if (*it == "qlbase")
        {
            const double ql_threshold = 0.;
            calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
            nerror += cross->cross_height_threshold(fields->atmp["tmp1"]->data, fields->atmp["tmp1"]->databot, fields->atmp["tmp2"]->data, grid->z, ql_threshold, Bottom_to_top, "qlbase");
        }
        else if (*it == "qltop")
        {
            const double ql_threshold = 0.;
            calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
            nerror += cross->cross_height_threshold(fields->atmp["tmp1"]->data, fields->atmp["tmp1"]->databot, fields->atmp["tmp2"]->data, grid->z, ql_threshold, Top_to_bottom, "qltop");
        }
        else if (*it == "maxthvcloud")
        {
            calc_liquid_water(fields->atmp["tmp1"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
            calc_maximum_thv_perturbation_cloud(fields->atmp["tmp2"]->databot, fields->atmp["tmp2"]->data,
                                                fields->sp["thl"]->data, fields->sp["qt"]->data, fields->atmp["tmp1"]->data, pref, fields->atmp["tmp2"]->datamean);
            nerror += cross->cross_plane(fields->atmp["tmp2"]->databot, fields->atmp["tmp1"]->data, "maxthvcloud");
        }
        else if (*it == "bbot" or *it == "bfluxbot")
        {
            calc_buoyancy_bot(fields->atmp["tmp1"]->data, fields->atmp["tmp1"]->databot, fields->sp[thvar ]->data, fields->sp[thvar]->databot, fields->sp["qt"]->data, fields->sp["qt"]->databot, thvref, thvrefh);
            calc_buoyancy_fluxbot(fields->atmp["tmp1"]->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);

            if (*it == "bbot")
                nerror += cross->cross_plane(fields->atmp["tmp1"]->databot, fields->atmp["tmp1"]->data, "bbot");
            else if (*it == "bfluxbot")
                nerror += cross->cross_plane(fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->data, "bfluxbot");
        }
        // BvS:micro 
        else if (*it == "qrpath")
        {
            nerror += cross->cross_path(fields->sp["qr"]->data, fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, "qrpath");
        }
    }

    if (nerror)
        throw 1;
}

void Thermo_moist::exec_dump()
{
    for (std::vector<std::string>::const_iterator it=dumplist.begin(); it<dumplist.end(); ++it)
    {
        // TODO BvS restore getThermoField(), the combination of checkThermoField with getThermoField is more elegant...
        if (*it == "b")
            calc_buoyancy(fields->atmp["tmp2"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, fields->atmp["tmp1"]->data, thvref);
        else if (*it == "ql")
            calc_liquid_water(fields->atmp["tmp2"]->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
        else
            throw 1;

        model->dump->save_dump(fields->atmp["tmp2"]->data, fields->atmp["tmp1"]->data, *it);
    }
}

bool Thermo_moist::check_field_exists(const std::string name)
{
    if (name == "b" || name == "ql")
        return true;
    else
        return false;
}

#ifndef USECUDA
void Thermo_moist::get_thermo_field(Field3d* fld, Field3d* tmp, const std::string name, bool cyclic)
{
    const int kcells = grid->kcells;

    // BvS: getThermoField() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
    // Pass dummy as rhoref,thvref to prevent overwriting base state
    double* restrict tmp2 = fields->atmp["tmp2"]->data;
    if (swupdatebasestate)
        calc_base_state(pref, prefh, &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], exnref, exnrefh,
                fields->sp[thvar]->datamean, fields->sp["qt"]->datamean);

    if (name == "b")
        calc_buoyancy(fld->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref, tmp->data, thvref);
    else if (name == "ql")
        calc_liquid_water(fld->data, fields->sp[thvar]->data, fields->sp["qt"]->data, pref);
    else if (name == "N2")
        calc_N2(fld->data, fields->sp[thvar]->data, grid->dzi, thvref);
    else
        throw 1;

    if (cyclic)
        grid->boundary_cyclic(fld->data);
}
#endif

#ifndef USECUDA
void Thermo_moist::get_buoyancy_surf(Field3d* bfield)
{
    calc_buoyancy_bot(bfield->data, bfield->databot,
                      fields->sp[thvar]->data, fields->sp[thvar]->databot,
                      fields->sp["qt"]->data, fields->sp["qt"]->databot,
                      thvref, thvrefh);
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

#ifndef USECUDA
void Thermo_moist::get_buoyancy_fluxbot(Field3d *bfield)
{
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp[thvar]->databot, fields->sp[thvar]->datafluxbot, fields->sp["qt"]->databot, fields->sp["qt"]->datafluxbot, thvrefh);
}
#endif

void Thermo_moist::get_prog_vars(std::vector<std::string> *list)
{
    list->push_back(thvar);
    list->push_back("qt");
}

double Thermo_moist::get_buoyancy_diffusivity()
{
    // Use the diffusivity from the liquid water potential temperature
    return fields->sp[thvar]->visc;
}

namespace
{
    inline double sat_adjust(const double thl, const double qt, const double p, const double exn)
    {
        int niter = 0;
        int nitermax = 30;
        double ql, tl, tnr_old = 1.e9, tnr, qs=0;
        tl = thl * exn;

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment
        if(qt-qsat(p, tl) <= 0)
            return 0.;

        tnr = tl;
        while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5 && niter < nitermax)
        {
            ++niter;
            tnr_old = tnr;
            qs = qsat(p,tnr);
            tnr = tnr - (tnr+(Lv/cp)*qs-tl-(Lv/cp)*qt)/(1+(std::pow(Lv,2)*qs)/ (Rv*cp*std::pow(tnr,2)));
        }

        if (niter == nitermax)
        {
            printf("Saturation adjustment not converged!! [thl=%f K, qt=%f kg/kg, p=%f p]\n",thl,qt,p);
            throw 1;
        }

        ql = std::max(0.,qt - qs);
        return ql;
    }
}

/**
 * This function calculates the hydrostatic pressure at full and half levels,
 * with option to return base state profiles like reference density and temperature
 * Solves: dpi/dz=-g/thv with pi=cp*(p/p0)**(rd/cp)
 * @param pref Pointer to output hydrostatic pressure array (full level)
 * @param prefh Pointer to output hydrostatic pressure array (half level)
 * @param rho Pointer to output density array (full level)
 * @param rhoh Pointer to output density array (half level)
 * @param thv Pointer to output virtual potential temperature array (full level)
 * @param thvh Pointer to output virtual potential temperature array (half level)
 * @param ex Pointer to output exner array (full level)
 * @param exh Pointer to output exner array (half level)
 * @param thlmean Pointer to input liq. water potential temperature array (horizontal mean, full level)
 * @param qtmean Pointer to input tot. moisture mix. ratio  array (horizontal mean, full level)
 */
void  Thermo_moist::calc_base_state(double* restrict pref,    double* restrict prefh,
                                    double* restrict rho,     double* restrict rhoh,
                                    double* restrict thv,     double* restrict thvh,
                                    double* restrict ex,      double* restrict exh,
                                    double* restrict thlmean, double* restrict qtmean)
{
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double thlsurf = interp2(thlmean[kstart-1], thlmean[kstart]);
    const double qtsurf  = interp2(qtmean[kstart-1],  qtmean[kstart]);

    double ql;

    // Calculate the values at the surface (half level == kstart)
    prefh[kstart] = pbot;
    exh[kstart]   = exner(prefh[kstart]);
    ql            = sat_adjust(thlsurf, qtsurf, prefh[kstart], exh[kstart]);
    thvh[kstart]  = virtual_temperature(exh[kstart], thlsurf, qtsurf, ql);
    rhoh[kstart]  = pbot / (Rd * exh[kstart] * thvh[kstart]);

    // Calculate the first full level pressure
    pref[kstart]  = prefh[kstart] * std::exp(-grav * grid->z[kstart] / (Rd * exh[kstart] * thvh[kstart]));

    for (int k=kstart+1; k<kend+1; ++k)
    {
        // 1. Calculate remaining values (thv and rho) at full-level[k-1]
        ex[k-1]  = exner(pref[k-1]);
        ql       = sat_adjust(thlmean[k-1], qtmean[k-1], pref[k-1], ex[k-1]);
        thv[k-1] = virtual_temperature(ex[k-1], thlmean[k-1], qtmean[k-1], ql);
        rho[k-1] = pref[k-1] / (Rd * ex[k-1] * thv[k-1]);

        // 2. Calculate pressure at half-level[k]
        prefh[k] = prefh[k-1] * std::exp(-grav * grid->dz[k-1] / (Rd * ex[k-1] * thv[k-1]));
        exh[k]   = exner(prefh[k]);

        // 3. Use interpolated conserved quantities to calculate half-level[k] values
        const double thli = interp2(thlmean[k-1], thlmean[k]);
        const double qti  = interp2(qtmean [k-1], qtmean [k]);
        const double qli  = sat_adjust(thli, qti, prefh[k], exh[k]);

        thvh[k]  = virtual_temperature(exh[k], thli, qti, qli);
        rhoh[k]  = prefh[k] / (Rd * exh[k] * thvh[k]);

        // 4. Calculate pressure at full-level[k]
        pref[k] = pref[k-1] * std::exp(-grav * grid->dzh[k] / (Rd * exh[k] * thvh[k]));
    }

    pref[kstart-1] = 2.*prefh[kstart] - pref[kstart];
}

void Thermo_moist::calc_buoyancy_tend_2nd(double* restrict wt, double* restrict thl, double* restrict qt,
                                          double* restrict ph, double* restrict thlh, double* restrict qth,
                                          double* restrict ql, double* restrict thvrefh)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double tl, exnh;

    for (int k=grid->kstart+1; k<grid->kend; k++)
    {
        exnh = exner(ph[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                thlh[ij] = interp2(thl[ijk-kk], thl[ijk]);
                qth[ij]  = interp2(qt[ijk-kk], qt[ijk]);
                tl       = thlh[ij] * exnh;
                // Calculate first estimate of ql using Tl
                // if ql(Tl)>0, saturation adjustment routine needed
                ql[ij]  = qth[ij]-qsat(ph[k],tl);
            }
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij  = i + j*jj;
                if (ql[ij]>0)   // already doesn't vectorize because of iteration in sat_adjust()
                {
                    ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh);
                }
                else
                    ql[ij] = 0.;
            }
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
            }
    }
}


void Thermo_moist::calc_buoyancy(double* restrict b, double* restrict thl, double* restrict qt,
                                 double* restrict p, double* restrict ql, double* restrict thvref)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double tl, ex;

    for (int k=0; k<grid->kcells; k++)
    {
        ex = exner(p[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                tl  = thl[ijk] * ex;
                ql[ij]  = qt[ijk]-qsat(p[k],tl);   // not real ql, just estimate
            }

        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                if (ql[ij] > 0)
                    ql[ij] = sat_adjust(thl[ijk], qt[ijk], p[k], ex);
                else
                    ql[ij] = 0.;
            }

        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;
                b[ijk] = buoyancy(ex, thl[ijk], qt[ijk], ql[ij], thvref[k]);
            }
    }

    grid->boundary_cyclic(b);
}

// Calculate the maximum in-cloud virtual temperature perturbation (thv - <thv>) with height
void Thermo_moist::calc_maximum_thv_perturbation_cloud(double* restrict thv_pert, double* restrict thv, double* restrict thl,
                                                       double* restrict qt, double* restrict ql, double* restrict p, double* restrict thvmean)
{
    double ex;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Set field to zero
    for (int j=grid->jstart; j<grid->jend; j++)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij = i + j*jj;
            thv_pert[ij] = NC_FILL_DOUBLE; // BvS not so nice..
        }

    // Calculate virtual temperature
    for (int k=grid->kstart; k<grid->kend; k++)
    {
        ex = exner(p[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                thv[ijk] = virtual_temperature(ex, thl[ijk], qt[ijk], ql[ijk]);
            }
    }

    // Calculate domain averaged virtual temperature
    grid->calc_mean(thvmean, thv, grid->kcells);

    // Calculate maximum in-cloud virtual temperature perturbation per column
    for (int k=grid->kstart; k<grid->kend; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + k*kk;

                if(ql[ijk] > 0)
                {
                    if(thv_pert[ij] == NC_FILL_DOUBLE)
                        thv_pert[ij] = thv[ijk] - thvmean[k];
                    else
                        thv_pert[ij] = std::max(thv_pert[ij], thv[ijk] - thvmean[k]);
                }
            }
}

void Thermo_moist::calc_liquid_water(double* restrict ql, double* restrict thl, double* restrict qt, double* restrict p)
{
    double ex;

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Fill ghost cells with zeros to prevent problems in calculating ql or qlcore masks
    for (int k=0; k<grid->kstart; k++)
    {
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ql[ijk] = 0.;
            }
    }

    for (int k=grid->kend; k<grid->kcells; k++)
    {
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ql[ijk] = 0.;
            }
    }

    // Calculate the ql field
    for (int k=grid->kstart; k<grid->kend; k++)
    {
        ex = exner(p[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ql[ijk] = sat_adjust(thl[ijk], qt[ijk], p[k], ex);
            }
    }
}

void Thermo_moist::calc_N2(double* restrict N2, double* restrict thl, double* restrict dzi,
                           double* restrict thvref)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                N2[ijk] = grav/thvref[k]*0.5*(thl[ijk+kk] - thl[ijk-kk])*dzi[k];
            }
}

void Thermo_moist::calc_buoyancy_bot(double* restrict b,      double* restrict bbot,
                                     double* restrict thl,    double* restrict thlbot,
                                     double* restrict qt,     double* restrict qtbot,
                                     double* restrict thvref, double* restrict thvrefh)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    // assume no liquid water at the lowest model level
    for (int j=0; j<grid->jcells; j++)
        #pragma ivdep
        for (int i=0; i<grid->icells; i++)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            bbot[ij ] = buoyancy_no_ql(thlbot[ij], qtbot[ij], thvrefh[kstart]);
            b   [ijk] = buoyancy_no_ql(thl[ijk], qt[ijk], thvref[kstart]);
        }
}

void Thermo_moist::calc_buoyancy_fluxbot(double* restrict bfluxbot, double* restrict thlbot, double* restrict thlfluxbot, double* restrict qtbot, double* restrict qtfluxbot,
        double* restrict thvrefh)
{
    const int jj = grid->icells;
    const int kstart = grid->kstart;

    // assume no liquid water at the lowest model level
    for (int j=0; j<grid->jcells; j++)
        #pragma ivdep
        for (int i=0; i<grid->icells; i++)
        {
            const int ij = i + j*jj;
            bfluxbot[ij] = buoyancy_flux_no_ql(thlbot[ij], thlfluxbot[ij], qtbot[ij], qtfluxbot[ij], thvrefh[kstart]);
        }
}

void Thermo_moist::init_stat()
{
    // Add variables to the statistics
    if (stats->get_switch() == "1")
    {
        /* Add fixed base-state density and temperature profiles. Density should probably be in fields (?), but
           there the statistics are initialized before thermo->create() is called */
        stats->add_fixed_prof("rhoref",  "Full level basic state density", "kg m-3", "z",  fields->rhoref );
        stats->add_fixed_prof("rhorefh", "Half level basic state density", "kg m-3", "zh", fields->rhorefh);
        stats->add_fixed_prof("thvref",  "Full level basic state virtual potential temperature", "K", "z", thvref );
        stats->add_fixed_prof("thvrefh", "Half level basic state virtual potential temperature", "K", "zh",thvrefh);

        if (swupdatebasestate == 1)
        {
            stats->add_prof("ph",   "Full level hydrostatic pressure", "Pa",     "z" );
            stats->add_prof("phh",  "Half level hydrostatic pressure", "Pa",     "zh");
            stats->add_prof("rho",  "Full level density",  "kg m-3", "z" );
            stats->add_prof("rhoh", "Half level density",  "kg m-3", "zh");
        }
        else
        {
            stats->add_fixed_prof("ph",  "Full level hydrostatic pressure", "Pa", "z",  pref );
            stats->add_fixed_prof("phh", "Half level hydrostatic pressure", "Pa", "zh", prefh);
        }

        stats->add_prof("b", "Buoyancy", "m s-2", "z");
        for (int n=2; n<5; ++n)
        {
            std::stringstream ss;
            ss << n;
            std::string sn = ss.str();
            stats->add_prof("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
        }

        stats->add_prof("bgrad", "Gradient of the buoyancy", "m s-3", "zh");
        stats->add_prof("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
        stats->add_prof("bdiff", "Diffusive flux of the buoyancy", "m2 s-3", "zh");
        stats->add_prof("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

        stats->add_prof("ql", "Liquid water mixing ratio", "kg kg-1", "z");
        stats->add_prof("cfrac", "Cloud fraction", "-","z");

        stats->add_time_series("lwp", "Liquid water path", "kg m-2");
        stats->add_time_series("ccover", "Projected cloud cover", "-");

        // BvS:micro 
        if(swmicro == "2mom_warm")
        {
            stats->add_time_series("rwp", "Rain water path", "kg m-2");

            if(swmicrobudget == "1")
            {
                stats->add_prof("auto_qrt" , "Autoconversion tendency qr", "kg kg-1 s-1", "z");
                stats->add_prof("auto_nrt" , "Autoconversion tendency nr", "m-3 s-1", "z");
                stats->add_prof("auto_thlt", "Autoconversion tendency thl", "K s-1", "z");
                stats->add_prof("auto_qtt" , "Autoconversion tendency qt", "kg kg-1 s-1", "z");

                stats->add_prof("evap_qrt" , "Evaporation tendency qr", "kg kg-1 s-1", "z");
                stats->add_prof("evap_nrt" , "Evaporation tendency nr", "m-3 s-1", "z");
                stats->add_prof("evap_thlt", "Evaporation tendency thl", "K s-1", "z");
                stats->add_prof("evap_qtt" , "Evaporation tendency qt", "kg kg-1 s-1", "z");

                stats->add_prof("accr_qrt" , "Accretion tendency qr", "kg kg-1 s-1", "z");
                stats->add_prof("accr_thlt", "Accretion tendency thl", "K s-1", "z");
                stats->add_prof("accr_qtt" , "Accretion tendency qt", "kg kg-1 s-1", "z");

                stats->add_prof("scbr_nrt" , "Selfcollection and breakup tendency nr", "m-3 s-1", "z");

                stats->add_prof("sed_qrt"  , "Sedimentation tendency qr", "kg kg-1 s-1", "z");
                stats->add_prof("sed_nrt"  , "Sedimentation tendency nr", "m-3 s-1", "z");
            }
        }
    }
}

void Thermo_moist::init_cross()
{
    if (model->cross->get_switch() == "1")
    {
        allowedcrossvars.push_back("b");
        allowedcrossvars.push_back("bbot");
        allowedcrossvars.push_back("bfluxbot");
        if (grid->swspatialorder == "4")
            allowedcrossvars.push_back("blngrad");
        allowedcrossvars.push_back("ql");
        allowedcrossvars.push_back("qlpath");
        allowedcrossvars.push_back("qlbase");
        allowedcrossvars.push_back("qltop");
        allowedcrossvars.push_back("maxthvcloud");  // yikes

        // BvS:micro 
        if(swmicro == "2mom_warm")
            allowedcrossvars.push_back("qrpath");

        // Get global cross-list from cross.cxx
        std::vector<std::string> *crosslist_global = model->cross->get_crosslist();

        // Check input list of cross variables (crosslist)
        std::vector<std::string>::iterator it=crosslist_global->begin();
        while (it != crosslist_global->end())
        {
            if (std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it))
            {
                // Remove variable from global list, put in local list
                crosslist.push_back(*it);
                crosslist_global->erase(it); // erase() returns iterator of next element..
            }
            else
                ++it;
        }

        // Sort crosslist to group ql and b variables
        std::sort(crosslist.begin(),crosslist.end());
    }
}

void Thermo_moist::init_dump()
{
    if (model->dump->get_switch() == "1")
    {
        // Get global cross-list from cross.cxx
        std::vector<std::string> *dumplist_global = model->dump->get_dumplist();

        // Check if fields in dumplist are retrievable thermo fields
        std::vector<std::string>::iterator dumpvar=dumplist_global->begin();
        while (dumpvar != dumplist_global->end())
        {
            if (check_field_exists(*dumpvar))
            {
                // Remove variable from global list, put in local list
                dumplist.push_back(*dumpvar);
                dumplist_global->erase(dumpvar); // erase() returns iterator of next element..
            }
            else
                ++dumpvar;
        }
    }
}

void Thermo_moist::calc_buoyancy_tend_4th(double* restrict wt, double* restrict thl,  double* restrict qt,
                                          double* restrict ph, double* restrict thlh, double* restrict qth,
                                          double* restrict ql, double* restrict thvrefh)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    double tl, exnh;

    for (int k=grid->kstart+1; k<grid->kend; k++)
    {
        exnh = exner(ph[k]);
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk1;
                const int ij  = i + j*jj;
                thlh[ij] = interp4(thl[ijk-kk2], thl[ijk-kk1], thl[ijk], thl[ijk+kk1]);
                qth[ij]  = interp4(qt[ijk-kk2],  qt[ijk-kk1],  qt[ijk],  qt[ijk+kk1]);
                tl       = thlh[ij] * exnh;
                // Calculate first estimate of ql using Tl
                // if ql(Tl)>0, saturation adjustment routine needed
                ql[ij]  = qth[ij]-qsat(ph[k],tl);
            }
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij = i + j*jj;
                if (ql[ij]>0)   // already doesn't vectorize because of iteration in sat_adjust()
                    ql[ij] = sat_adjust(thlh[ij], qth[ij], ph[k], exnh);
                else
                    ql[ij] = 0.;
            }
        for (int j=grid->jstart; j<grid->jend; j++)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk1;
                const int ij  = i + j*jj;
                wt[ijk] += buoyancy(exnh, thlh[ij], qth[ij], ql[ij], thvrefh[k]);
            }
    }
}

