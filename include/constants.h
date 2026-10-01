/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants
{
    template<typename TF> constexpr TF kappa = 0.4;           // von Karman constant
    template<typename TF> constexpr TF grav  = 9.81;          // Gravitational acceleration [m s-2]
    template<typename TF> constexpr TF e_rot = 7.2921e-5;     // Earth rotation rate [s-1]
    template<typename TF> constexpr TF Rd    = 287.04;        // Gas constant for dry air [J K-1 kg-1]
    template<typename TF> constexpr TF Rv    = 461.5;         // Gas constant for water vapor [J K-1 kg-1]
    template<typename TF> constexpr TF cp    = 1005;          // Specific heat of air at constant pressure [J kg-1 K-1]
    template<typename TF> constexpr TF Lv    = 2.501e6;       // Latent heat of vaporization [J kg-1]
    template<typename TF> constexpr TF Lf    = 3.337e5;       // Latent heat of fusion [J kg-1]
    template<typename TF> constexpr TF Ls    = Lv<TF>+Lf<TF>; // Latent heat of sublimation [J kg-1]
    template<typename TF> constexpr TF T0    = 273.15;        // Freezing / melting temperature [K]
    template<typename TF> constexpr TF p0    = 1.e5;          // Reference pressure [Pa]
    template<typename TF> constexpr TF ep    = Rd<TF>/Rv<TF>;
    template<typename TF> constexpr TF rho_w = 1.e3;          // Density of water [kg m-3]
    template<typename TF> constexpr TF rho_i = 7.e2;          // Density of ice   [kg m-3]
    template<typename TF> constexpr TF mu0_min = 1e-6;        // Minimum value used for cos(sza)
    template<typename TF> constexpr TF sigma_b = 5.67e-8;     // Boltzmann constant [W m-1 K-1]
    template<typename TF> constexpr TF xmair = 28.9647;       // Molar mass of dry air [kg kmol-1]
    template<typename TF> constexpr TF xmh2o = 18.01528;      // Molar mass of h2o [kg kmol-1]

    // Soil / land-surface specific constants
    template<typename TF> constexpr TF rho_C_matrix   = 1.6e6;   // Volumetric soil heat capacity [J m-3 K-1]
    template<typename TF> constexpr TF rho_C_water    = 4.18e6;  // Volumetric water heat capacity [J m-3 K-1]
    template<typename TF> constexpr TF gamma_T_matrix = 3.4293695508945325; //  pow(7.7, 0.4)*pow(2,(1-0.4));  // Heat conductivity soil [J s-1 m-1 K-1]
    template<typename TF> constexpr TF gamma_T_water  = 0.57;    // Heat conductivity water [J s-1 m-1 K-1]
    template<typename TF> constexpr TF wlmax = 0.0002;           // Max water depth on vegetation/soil

    // Limits on Obukhov length:
    template<typename TF> constexpr TF zL_max = 10.;
    template<typename TF> constexpr TF zL_min = -1.e4;

    // NOTE: we need to keep `sgstke_min` around epsilon<float>, which is ~1e-7.
    template<typename TF> constexpr TF sgstke_min = 1e-7;   // Minimum value SGS TKE [m2 s-2]

    // BvS: Perhaps better off back in defines.h? Separate namespace?
    const double        dtiny  = 1.e-30;
    const double        dsmall = 1.e-9;
    const double        dbig   = 1.e9;
    const double        dhuge  = 1.e30;
    const unsigned long ulhuge = ~0UL; //== ULONG_MAX;
}
#endif
