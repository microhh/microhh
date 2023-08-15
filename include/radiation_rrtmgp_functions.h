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

#ifndef RADIATION_RRTMGP_FUNCTIONS_H
#define RADIATION_RRTMGP_FUNCTIONS_H

//// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__
#else
#  define CUDA_MACRO
#endif


#include <iostream>
#include <iomanip>

#include "types.h"

namespace Radiation_rrtmgp_functions
{
    inline std::pair<Float, Float> calc_cos_zenith_angle(
            const Float lat, const Float lon, const int day_of_year,
            const Float seconds_since_midnight, const int year)
    {
        /* Based on: Paltridge, G. W. and Platt, C. M. R. (1976).
                     Radiative Processes in Meteorology and Climatology.
                     Elsevier, New York, 318 pp. */

        const Float pi = Float(M_PI);

        // Account for leap year
        int days_per_year;
        if ((year%4 == 0) && ((year%100 != 0) || (year%400 == 0)))
            days_per_year = 366;
        else
            days_per_year = 365;

        // DOY in time calculations are zero based:
        const int doy = day_of_year-1;

        // Lat/lon in radians
        const Float radlat = lat * pi/Float(180);
        const Float radlon = lon * pi/Float(180);

        // DOY in range (0,2*pi)
        const Float doy_pi = Float(2)*pi*doy/days_per_year;

        // Solar declination angle
        const Float declination_angle = \
                Float(0.006918) - Float(0.399912) * std::cos(doy_pi) + Float(0.070257) * std::sin(doy_pi)
              - Float(0.006758) * std::cos(Float(2)*doy_pi) + Float(0.000907) * std::sin(Float(2)*doy_pi)
              - Float(0.002697) * std::cos(Float(3)*doy_pi) + Float(0.00148)  * std::sin(Float(3)*doy_pi);

        // Hour angle in radians, using true solar time
        const Float a1 = (Float(1.00554) * doy - Float( 6.28306)) * pi/Float(180);
        const Float a2 = (Float(1.93946) * doy + Float(23.35089)) * pi/Float(180);
        const Float a3 = (Float(7.67825) * std::sin(a1) + Float(10.09176) * std::sin(a2)) / Float(60);

        const Float hour_solar_time = (seconds_since_midnight/Float(3600)) - a3 + radlon * (Float(180.)/pi/Float(15.0));
        const Float hour_angle = (hour_solar_time-Float(12))*Float(15.0)*(pi/Float(180));

        // Cosine of solar zenith angle
        const Float cos_zenith = std::sin(radlat) * std::sin(declination_angle)
                               + std::cos(radlat) * std::cos(declination_angle) * std::cos(hour_angle);

        const Float cos_elevation = std::cos(Float(0.5)*pi - std::acos(cos_zenith));

        // put maximum at -1 to prevent problems in single precision
        const Float cos_azimuth = std::max(Float(-1.), (
              std::cos(radlat) * std::sin(declination_angle)
            - std::sin(radlat) * std::cos(declination_angle) * std::cos(hour_angle) ) / cos_elevation);

        const Float azimuth = (hour_angle <= Float(0.)) ? std::acos(cos_azimuth) : Float(2.)*pi - std::acos(cos_azimuth);

        return std::pair<Float, Float>(cos_zenith, azimuth);
    }


    inline Float calc_sun_distance_factor(const Float frac_doy)
    {
        // Based on: An Introduction to Atmospheric Radiation, Liou, Eq. 2.2.9.
        constexpr Float an [] = {1.000110, 0.034221, 0.000719};
        constexpr Float bn [] = {0,        0.001280, 0.000077};

        const Float pi = Float(M_PI);

        // NOTE: DOY in time calculations are zero based
        const Float t = Float(2)*pi*(frac_doy-Float(1))/Float(365);

        Float factor = Float(0);
        for (int n=0; n<3; ++n)
            factor += an[n]*std::cos(Float(n)*t) + bn[n]*std::sin(Float(n)*t);

        return factor;
    }


}
#endif
