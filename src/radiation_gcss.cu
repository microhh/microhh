/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 * Copyright (c) 2018-2019 Elynn Wu
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

#include "tools.h"
#include "radiation_gcss.h"

namespace
{
     template<typename TF>
     TF calc_zenith(struct tm datetime, TF lat, TF lon)
     {
         const TF pi        = M_PI;
         const TF year2days = 365.;
         const TF piAngle   = 180.;
         const TF day2secs  = 86400.;
         const TF z1        = 279.934;
         const TF z2        = 1.914827;
         const TF z3        = 0.7952;
         const TF z4        = 0.019938;
         const TF z5        = 0.00162;
         const TF z6        = 23.4439;
         TF time2sec = (datetime.tm_yday + 1) +
                            lon / 360. +
                           (datetime.tm_hour * 3600. +
                            datetime.tm_min * 60. +
                            datetime.tm_sec) / day2secs;
         TF day    = floor(time2sec);
         TF lamda  = lat * pi / piAngle;
         TF d      = 2. * pi * int(time2sec) / year2days;
         TF sig    = d + pi/piAngle * (z1 + z2*std::sin(d)
                                                   - z3*std::cos(d)
                                                   + z4*std::sin(2.*d)
                                                   - z5*std::cos(2.*d));
         TF del     = std::asin(std::sin(z6*pi / piAngle)*std::sin(sig));
         TF h       = 2. * pi * ((time2sec - day) - 0.5);
         TF mu      = std::sin(lamda) * std::sin(del) + std::cos(lamda) * std::cos(del) * std::cos(h);
         return mu;
     }

     template<typename TF> __device__
     void sunray(const TF mu, const int i, const int j,
         const int kstart, const int kend, const int icells, const int ijcells,
         TF* __restrict__  tau, const TF tauc,
         TF* __restrict__  swn)
     {
         const int jj = icells;
         const int kk = ijcells;
         TF o_c1 = TF(0.9);
         TF o_c2 = TF(2.75);
         TF o_c3 = TF(0.09);
         TF sw0 = TF(1100.);
         TF gc  = TF(0.85);
         TF sfc_albedo = TF(0.05);
         TF taucde = TF(0.);
         TF taupath = TF(0.);
         TF omega  = TF(1.) - TF(1.e-3) * (o_c1 + o_c2 * (mu+TF(1.)) * exp(-o_c3 * tauc)); //fouquart
         TF ff     = gc * gc;
         TF gcde   = gc / (TF(1.) + gc);
         taucde = ( TF(1.0) - omega*ff) * tauc;
         TF omegade = (TF(1.)-ff) * omega/(TF(1.) - omega*ff);
         TF x1  = TF(1.) - omegade * gcde;
         TF x2  = TF(1.) - omegade;
         TF rk  = sqrt(TF(3.) * x2 * x1);
         TF mu2 = mu * mu;
         TF x3  = TF(4.) * (TF(1.) - rk*rk*mu2);
         TF rp  = sqrt(TF(3.) * x2/x1);
         TF alpha = TF(3.) * omegade * mu2 * (TF(1.) + gcde*x2) / x3;
         TF beta  = TF(3.) * omegade * mu * (TF(1.) + TF(3.)*gcde*mu2*x2) / x3;

         TF rtt = TF(2.0/3.0);
         TF exmu0 = exp(-taucde / mu);
         TF expk  = exp(rk * taucde);
         TF exmk  = TF(1.) / expk;
         TF xp23p = TF(1.) + rtt*rp;
         TF xm23p = TF(1.) - rtt*rp;
         TF ap23b = alpha + rtt*beta;

         TF t1 = TF(1.) - sfc_albedo - rtt * (TF(1.) + sfc_albedo) * rp;
         TF t2 = TF(1.) - sfc_albedo + rtt * (TF(1.) + sfc_albedo) * rp;
         TF t3 = (TF(1.) - sfc_albedo) * alpha - rtt * (TF(1.) + sfc_albedo) * beta + sfc_albedo*mu;
         TF c2 = (xp23p*t3*exmu0 - t1*ap23b*exmk) / (xp23p*t2*expk - xm23p*t1*exmk);
         TF c1 = (ap23b - c2*xm23p)/xp23p;


         for (int k=kend-1;k>=kstart;--k)
         {
             const int ijk  = i + j*jj + k*kk;
             taupath += ( TF(1.) - omega*ff ) * tau[ijk];
             swn[ijk] = sw0 * TF(4./3.) * (rp * (c1*exp(-rk*taupath)
             - c2 * exp(rk*taupath)) - beta * exp(-taupath/mu))
             + mu * sw0 * exp(-taupath / mu);
         }
     }

     template<typename TF> __global__
     void calc_gcss_rad_SW_g(TF* __restrict__ swn, TF* __restrict__ ql,TF* __restrict__ qt,
                   TF* __restrict__ tau,
                   TF* __restrict__ rhoref, TF mu,
                   TF* __restrict__ z, TF* __restrict__ dz,
                   const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                   const int icells, const int ijcells, const int ncells)
     {
         const int jj = icells;
         const int kk = ijcells;
         const TF rho_l = 1000.;
         const TF reff = 1.E-5;
         const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
         const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

         if (i < iend && j < jend)
         {
             TF tauc = 0.;
             for (int k=kstart;k<kend;++k)
             {
                 const int ijk  = i + j*jj + k*kk;
                 tau[ijk] = TF(1.5) * ql[ijk] * rhoref[k] * dz[k] / reff / rho_l;
                 tauc = tauc + tau[ijk];
             }
             sunray<TF>(mu, i, j,
                 kstart, kend, icells, ijcells,
                 tau, tauc, swn);
         }
     }
     template<typename TF> __global__
     void calc_gcss_rad_LW_g(TF* __restrict__ flx, TF* __restrict__ ql,TF* __restrict__ qt, TF* __restrict__ lwp,
                   TF* __restrict__ rhoref, TF fr0, TF fr1, TF xka, TF div, TF cp,
                   TF* __restrict__ z, TF* __restrict__ dz,
                   const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                   const int jj, const int kk)
     {
         const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
         const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

         if (i < iend && j < jend)
         {
             const int ij  = i + j*jj;
             int ki = kend; //pblh index
             TF fact;
             lwp[ij] = TF(0.0);

             ki = kend; //set to top of domain
             for (int k=kstart; k<kend; ++k)
             {
                 const int ijk  = i + j*jj + k*kk;

                 lwp[ij] += ql[ijk] *rhoref[k] * dz[k];
                 flx[ijk] = fr1 * exp(TF(-1.0) * xka * lwp[ij]);
                 if ( (ql[ijk] > TF(0.01E-3) ) && ( qt[ijk] >= TF(0.008) ) ) ki = k; //this is the PBLH index
             }
             fact = max(TF(0.), div * cp * rhoref[ki]);

             flx[ij + kstart*kk] += fr0 * exp(TF(-1.0) * xka *lwp[ij]);
             for (int k=kstart+1; k<kend; ++k)
             {
                 const int ijk  = i + j*jj + k*kk;

                 lwp[ij]  -= ql[ijk] *rhoref[k] * dz[k];
                 flx[ijk] += fr0 * exp(TF(-1.0) * xka * lwp[ij]);
             }


             for (int k=ki; k<kend; ++k)
             {
                 const int ijk  = i + j*jj + k*kk;
                 flx[ijk] += fact * ( TF(0.25) * pow(z[k]-z[ki],TF(1.333)) + z[ki] * pow(z[k]-z[ki],TF(0.33333)) );
             }

         }
     }

     template<typename TF> __global__
     void update_temperature(TF* __restrict__ tt, TF* __restrict__ flx, TF cp, TF* __restrict__ rhoref,
                             TF* dzi, int istart, int jstart, int kstart,
                             int iend,   int jend,   int kend, int jj, int kk)
     {
         const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
         const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
         const int k = blockIdx.z + kstart + 1;

         if (i < iend && j < jend && k < kend)
         {
             const int ijk = i + j*jj + k*kk;
             tt[ijk] += - (flx[ijk+kk] - flx[ijk]) * dzi[k] / (rhoref[k] * cp);
         }
     }
    }

    #ifdef USECUDA
    template<typename TF>
    void Radiation_gcss<TF>::exec(Thermo<TF>& thermo, double time, Timeloop<TF>& timeloop)
    {
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;

    int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    auto tmp = fields.get_tmp_g();
    auto flx = fields.get_tmp_g();
    auto ql  = fields.get_tmp_g();

    thermo.get_thermo_field_g(*ql,"ql",false);

    calc_gcss_rad_LW_g<<<gridGPU, blockGPU>>>(flx->fld_g, ql->fld_g, fields.sp.at("qt")->fld_g,
        tmp->fld_g,fields.rhoref_g, fr0, fr1, xka, div, Constants::cp<TF>,gd.z_g, gd.dz_g,
        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
        gd.icells, gd.ijcells);

    update_temperature<<<gridGPU, blockGPU>>>(
            fields.st.at("thl")->fld_g, flx->fld_g, Constants::cp<TF>, fields.rhoref_g,
            gd.dzi_g, gd.istart,  gd.jstart, gd.kstart,
            gd.iend,  gd.jend,   gd.kend, gd.icells, gd.ijcells);

    struct tm current_datetime;
    current_datetime = timeloop.get_phytime();
    TF mu = calc_zenith(current_datetime, lat, lon);
    if (mu > mu_min)
    {
        calc_gcss_rad_SW_g<<<gridGPU, blockGPU>>>(flx->fld_g, ql->fld_g, fields.sp.at("qt")->fld_g,
            tmp->fld_g, fields.rhoref_g, mu, gd.z_g, gd.dzi_g,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells, gd.ncells);

        update_temperature<<<gridGPU, blockGPU>>>(
            fields.st.at("thl")->fld_g, flx->fld_g, Constants::cp<TF>, fields.rhoref_g,
            gd.dzi_g, gd.istart,  gd.jstart, gd.kstart,
            gd.iend,  gd.jend,   gd.kend, gd.icells, gd.ijcells);
    }
    fields.release_tmp_g(tmp);
    fields.release_tmp_g(flx);
    fields.release_tmp_g(ql);
}
#endif

#ifdef USECUDA
template<typename TF>
void Radiation_gcss<TF>::exec_column(Column<TF>& column, Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    const TF no_offset = 0.;
    auto flx = fields.get_tmp_g();

    get_radiation_field_g(*flx,"lflx",thermo, timeloop);
    column.calc_column("lflx", flx->fld_g, no_offset);
    get_radiation_field_g(*flx,"sflx",thermo,timeloop);
    column.calc_column("sflx", flx->fld_g, no_offset);

    fields.release_tmp_g(flx);
}
#endif


#ifdef USECUDA
template<typename TF>
void Radiation_gcss<TF>::get_radiation_field_g(Field3d<TF>& fld, std::string name, Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;

    int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    if (name == "lflx")
    {
        auto lwp = fields.get_tmp_g();
        auto ql  = fields.get_tmp_g();
        thermo.get_thermo_field_g(*ql,"ql",false);

        calc_gcss_rad_LW_g<<<gridGPU, blockGPU>>>(fld.fld_g, ql->fld_g,fields.sp.at("qt")->fld_g, lwp->fld_g,
            fields.rhoref_g, fr0, fr1, xka, div, Constants::cp<TF>,gd.z_g, gd.dz_g,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
        fields.release_tmp_g(lwp);
        fields.release_tmp_g(ql);
    }

    else if (name == "sflx")
    {
        struct tm current_datetime;
        current_datetime = timeloop.get_phytime();
        TF mu = calc_zenith(current_datetime, lat, lon);
        if (mu > mu_min) //if daytime, call SW (make a function for day/night determination)
        {
            auto ql  = fields.get_tmp_g();
            auto tau = fields.get_tmp_g();
            thermo.get_thermo_field_g(*ql,"ql",false);
            calc_gcss_rad_SW_g<<<gridGPU, blockGPU>>>(fld.fld_g, ql->fld_g, fields.sp.at("qt")->fld_g,
                tau->fld_g, fields.rhoref_g, mu, gd.z_g, gd.dzi_g,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells, gd.ncells);

            fields.release_tmp_g(ql);
            fields.release_tmp_g(tau);
        }
        else
        {
            const int nblock = 256;
            const int ngrid  = gd.ncells/nblock + (gd.ncells%nblock > 0);
            set_to_val<<<ngrid, nblock>>>(fld.fld_g, gd.ncells, TF(0.));
        }
    }
    else
    {
        std::string msg = "get_radiation_field \"" + name + "\" not supported";
        throw std::runtime_error(msg);
    }
}
#endif
template class Radiation_gcss<double>;
template class Radiation_gcss<float>;
