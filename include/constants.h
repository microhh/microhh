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

#ifndef CONSTANTS
#define CONSTANTS

#include <climits>

namespace Constants
{
    const double kappa = 0.4;        ///< von Karman constant
    const double grav  = 9.81;       ///< Gravitational acceleration [m s-2]
    const double Rd    = 287.04;     ///< Gas constant for dry air [J K-1 kg-1] 
    const double Rv    = 461.5;      ///< Gas constant for water vapor [J K-1 kg-1]
    const double cp    = 1005;       ///< Specific heat of air at constant pressure [J kg-1 K-1]
    const double Lv    = 2.5e6;      ///< Latent heat of condensation or vaporization [J kg-1]
    const double T0    = 273.15;     ///< Freezing / melting temperature [K]
    const double p0    = 1.e5;       ///< Reference pressure [pa]

    const double ep    = Rd/Rv;

    // Coefficients saturation vapor pressure estimation
    // Original MicroHH (/ UCLA-LES)
    const double c0 = 0.6105851e+03; 
    const double c1 = 0.4440316e+02; 
    const double c2 = 0.1430341e+01; 
    const double c3 = 0.2641412e-01; 
    const double c4 = 0.2995057e-03; 
    const double c5 = 0.2031998e-05; 
    const double c6 = 0.6936113e-08; 
    const double c7 = 0.2564861e-11; 
    const double c8 = -.3704404e-13; 

    // Coefficients Taylor expansion Arden Buck equation (1981) around T=T0
    const double c00  = +6.1121000000E+02;
    const double c10  = +4.4393067270E+01;
    const double c20  = +1.4279398448E+00;
    const double c30  = +2.6415206946E-02;
    const double c40  = +3.0291749160E-04;
    const double c50  = +2.1159987257E-06;
    const double c60  = +7.5015702516E-09;
    const double c70  = -1.5604873363E-12;
    const double c80  = -9.9726710231E-14;
    const double c90  = -4.8165754883E-17;
    const double c100 = +1.3839187032E-18;

    // Coefficients exner function estimation
    const double ex1 = 2.85611940298507510698e-06;
    const double ex2 = -1.02018879928714644313e-11;
    const double ex3 = 5.82999832046362073082e-17;
    const double ex4 = -3.95621945728655163954e-22;
    const double ex5 = 2.93898686274077761686e-27;
    const double ex6 = -2.30925409555411170635e-32;
    const double ex7 = 1.88513914720731231360e-37;

    // BvS: Perhaps better off back in defines.h? Separate namespace?
    const double        dtiny  = 1.e-30;
    const double        dsmall = 1.e-9;
    const double        dbig   = 1.e9;
    const double        dhuge  = 1.e30;
    const unsigned long ulhuge = ULONG_MAX;
}
#endif
