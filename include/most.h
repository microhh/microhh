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

#ifndef MOST

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

namespace most
{
  //
  // GRADIENT FUNCTIONS
  //
  CUDA_MACRO inline double phim_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return std::pow(1. + 3.6*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }

  CUDA_MACRO inline double phim_stable(const double zeta)
  {
    // Hogstrom, 1988
    return 1. + 4.8*zeta;
  }

  CUDA_MACRO inline double phim(const double zeta)
  {
    return (zeta <= 0.) ? phim_unstable(zeta) : phim_stable(zeta);
  }

  CUDA_MACRO inline double phih_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return std::pow(1. + 7.9*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }

  CUDA_MACRO inline double phih_stable(const double zeta)
  {
    // Hogstrom, 1988
    return 1. + 7.8*zeta;
  }

  CUDA_MACRO inline double phih(const double zeta)
  {
    return (zeta <= 0.) ? phih_unstable(zeta) : phih_stable(zeta);
  }

  //
  // INTEGRATED FUNCTIONS
  //
  CUDA_MACRO inline double psim_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return 3.*std::log( ( 1. + 1./phim_unstable(zeta) ) / 2.);
  }

  CUDA_MACRO inline double psim_stable(const double zeta)
  {
    // Hogstrom, 1988
    return -4.8*zeta;
  }

  CUDA_MACRO inline double psih_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return 3. * std::log( ( 1. + 1. / phih_unstable(zeta) ) / 2.);
  }

  CUDA_MACRO inline double psih_stable(const double zeta)
  {
    // Hogstrom, 1988
    return -7.8*zeta;
  }

  CUDA_MACRO inline double fm(const double zsl, const double z0m, const double L)
  {
    return (L <= 0.)
      ? constants::kappa / (std::log(zsl/z0m) - psim_unstable(zsl/L) + psim_unstable(z0m/L))
      : constants::kappa / (std::log(zsl/z0m) - psim_stable  (zsl/L) + psim_stable  (z0m/L));
  }

  CUDA_MACRO inline double fh(const double zsl, const double z0h, const double L)
  {
    return (L <= 0.)
      ? constants::kappa / (std::log(zsl/z0h) - psih_unstable(zsl/L) + psih_unstable(z0h/L))
      : constants::kappa / (std::log(zsl/z0h) - psih_stable  (zsl/L) + psih_stable  (z0h/L));
  }
}
#endif
