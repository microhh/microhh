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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "boundary_cyclic.h"
#include "thermo_moist_functions.h"
#include "constants.h"
#include "tools.h"

#include "microphys.h"
#include "microphys_2mom_warm.h"

using namespace Constants;
using namespace Thermo_moist_functions;
using namespace Micro_2mom_warm_constants;
using namespace Micro_2mom_warm_functions;

namespace
{
    template<typename TF> __global__
    void remove_negative_values_g(TF* __restrict__ field,
                                  int istart, int jstart, int kstart,
                                  int iend,   int jend,   int kend,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            field[ijk] = fmax(field[ijk], TF(0));
        }
    }

    // Autoconversion: formation of rain drop by coagulating cloud droplets
    template<typename TF> __global__
    void autoconversion_g(TF* const __restrict__ qrt, TF* const __restrict__ nrt,
                          TF* const __restrict__ qtt, TF* const __restrict__ thlt,
                          const TF* const __restrict__ qr,  const TF* const __restrict__ ql,
                          const TF* const __restrict__ rho, const TF* const __restrict__ exner,
                          const int istart, const int jstart, const int kstart,
                          const int iend,   const int jend,   const int kend,
                          const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        // BvS with the cpu and cuda implementation, these should really go to the header
        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES
        const TF k_cc   = 9.44e9;        // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48
        const TF nu_c   = 1;             // SB06, Table 1., same as UCLA-LES
        const TF kccxs  = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / pow(nu_c+1, TF(2));

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            if (ql[ijk] > ql_min<TF>)
            {
                const TF xc      = rho[k] * ql[ijk] / Nc0<TF>;    // Mean mass of cloud drops [kg]
                const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall);    // SB06, Eq 5
                const TF phi_au  = TF(600.) * pow(tau, TF(0.68)) * pow(TF(1.) - pow(tau, TF(0.68)), TF(3));    // UCLA-LES
                //const TF phi_au  = 400. * pow(tau, 0.7) * pow(1. - pow(tau, 0.7), 3);    // SB06, Eq 6
                const TF au_tend = rho[k] * kccxs * pow(ql[ijk], TF(2)) * pow(xc, TF(2)) *
                                       (TF(1.) + phi_au / pow(TF(1.)-tau, TF(2))); // SB06, eq 4

                qrt[ijk]  += au_tend;
                nrt[ijk]  += au_tend * rho[k] / x_star;
                qtt[ijk]  -= au_tend;
                thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * au_tend;
            }
        }
    }

    // Accreation: growth of raindrops collecting cloud droplets
    template<typename TF> __global__
    void accretion_g(TF* const __restrict__ qrt, TF* const __restrict__ qtt, TF* const __restrict__ thlt,
                     const TF* const __restrict__ qr,  const TF* const __restrict__ ql,
                     const TF* const __restrict__ rho, const TF* const __restrict__ exner,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF k_cr  = 5.25; // SB06, p49

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            if (ql[ijk] > ql_min<TF> && qr[ijk] > qr_min<TF>)
            {
                const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk]); // SB06, Eq 5
                const TF phi_ac  = pow(tau / (tau + TF(5e-5)), TF(4)); // SB06, Eq 8
                const TF ac_tend = k_cr * ql[ijk] *  qr[ijk] * phi_ac * pow(rho_0<TF> / rho[k], TF(0.5)); // SB06, Eq 7

                qrt[ijk]  += ac_tend;
                qtt[ijk]  -= ac_tend;
                thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * ac_tend;
            }
        }
    }




}

#ifdef USECUDA
template<typename TF>
void Microphys_2mom_warm<TF>::exec(Thermo<TF>& thermo, const double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    // Remove negative values from the qr and nr fields
    remove_negative_values_g<<<gridGPU, blockGPU>>>(fields.ap.at("qr")->fld_g,
        gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    cuda_check_error();

    remove_negative_values_g<<<gridGPU, blockGPU>>>(fields.ap.at("nr")->fld_g,
        gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    cuda_check_error();

    // Get cloud liquid water from thermodynamics
    auto ql = fields.get_tmp_g();
    thermo.get_thermo_field_g(*ql, "ql", false);

    // Get GPU pointers basestate pressure and exner from thermo
    TF* p = thermo.get_basestate_fld_g("pref");
    TF* exner = thermo.get_basestate_fld_g("exner");

    // ---------------------------------
    // Calculate microphysics tendencies
    // ---------------------------------
    
    // Autoconversion; formation of rain drop by coagulating cloud droplets
    autoconversion_g<<<gridGPU, blockGPU>>>(
        fields.st.at("qr")->fld_g, fields.st.at("nr")->fld_g,
        fields.st.at("qt")->fld_g, fields.st.at("thl")->fld_g,
        fields.sp.at("qr")->fld_g, ql->fld_g, fields.rhoref_g, exner,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);

    // Accretion; growth of raindrops collecting cloud droplets
    accretion_g<<<gridGPU, blockGPU>>>(
        fields.st.at("qr")->fld_g, fields.st.at("qt")->fld_g, fields.st.at("thl")->fld_g,
        fields.sp.at("qr")->fld_g, ql->fld_g, fields.rhoref_g, exner,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);

    fields.release_tmp_g(ql);
}
#endif

#ifdef USECUDA
template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}
#endif

template class Microphys_2mom_warm<double>;
template class Microphys_2mom_warm<float>;
