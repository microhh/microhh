#
#  MicroHH
#  Copyright (c) 2011-2024 Chiel van Heerwaarden
#  Copyright (c) 2011-2024 Thijs Heus
#  Copyright (c) 2014-2024 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

# Standard library

# Third-party.
import numpy as np
from numba import njit

# Local library
import puhhpy.constants as cst


@njit
def pow2(x):
    return x*x

@njit
def exner(p):
    return (p/cst.p0)**(cst.Rd/cst.cp)

@njit
def esat_liq(T):
    x = np.maximum(-75., T-cst.T0);
    return 611.21*np.exp(17.502*x / (240.97+x))

@njit
def esat_ice(T):
    x = np.maximum(-100., T-cst.T0)
    return 611.15*np.exp(22.452*x / (272.55+x))

@njit
def qsat_liq(p, T):
    return cst.ep*esat_liq(T)/(p-(1.-cst.ep)*esat_liq(T))

@njit
def qsat_ice(p, T):
    return cst.ep*esat_ice(T)/(p-(1.-cst.ep)*esat_ice(T))

@njit
def dqsatdT_liq(p, T):
    den = p - esat_liq(T)*(1.-cst.ep)
    return (cst.ep/den - (1. + cst.ep)*cst.ep*esat_liq(T)/pow2(den)) * cst.Lv*esat_liq(T) / (cst.Rv*pow2(T))

@njit
def dqsatdT_ice(p, T):
    den = p - esat_ice(T)*(1.-cst.ep)
    return (cst.ep/den + (1. - cst.ep)*cst.ep*esat_ice(T)/pow2(den)) * cst.Ls*esat_ice(T) / (cst.Rv*pow2(T))

@njit
def water_fraction(T):
    return max(0., min((T - 233.15) / (cst.T0 - 233.15), 1.))

@njit
def qsat(p, T):
    alpha = water_fraction(T)
    return alpha*qsat_liq(p, T) + (1.-alpha)*qsat_ice(p, T)

@njit
def virtual_temperature(exn, thl, qt, ql, qi):
    th = thl + cst.Lv*ql/(cst.cp*exn) + cst.Ls*qi/(cst.cp*exn)
    return th * (1. - (1. - cst.Rv/cst.Rd)*qt - cst.Rv/cst.Rd*(ql+qi));

@njit
def sat_adjust(thl, qt, p, use_ice=True):
    """
    Saturation adjustment of input `thl` and `qt` given `p`.
    The output T, ql, qi, qsat are such that qt = qsat(T,p) + ql + qi.
    """

    niter = 0
    nitermax = 10
    tnr_old = 1e9

    exn = exner(p)
    tl = thl*exn
    qs = qsat(p, tl)

    if (qt-qs <= 0.):
        return tl, 0., 0., qs

    # Warm adjustment.
    tnr = tl
    if tnr >= cst.T0 or not use_ice:
        while abs(tnr - tnr_old) / tnr_old > 1.e-5 and niter < nitermax:
            niter += 1
            tnr_old = tnr
            qs = qsat_liq(p, tnr)
            f = tnr - tl - cst.Lv/cst.cp*(qt - qs)
            f_prime = 1. + cst.Lv/cst.cp*dqsatdT_liq(p, tnr)
            tnr -= f / f_prime

        qs = qsat_liq(p, tnr)
        ql = max(0., qt - qs)
        return tnr, ql, 0., qs

    # Cold adjustment.
    else:
        while abs(tnr - tnr_old) / tnr_old > 1.e-5 and niter < nitermax:
            niter += 1
            tnr_old = tnr
            qs = qsat(p, tnr)

            alpha_w = water_fraction(tnr)
            alpha_i = 1. - alpha_w
            dalphadT = 0.025 if alpha_w > 0. and alpha_w < 1. else 0.

            dqsatdT_w = dqsatdT_liq(p, tnr)
            dqsatdT_i = dqsatdT_ice(p, tnr)

            f = tnr - tl - alpha_w*cst.Lv/cst.cp*qt - alpha_i*cst.Ls/cst.cp*qt \
                         + alpha_w*cst.Lv/cst.cp*qs + alpha_i*cst.Ls/cst.cp*qs

            f_prime = 1. \
                - dalphadT*cst.Lv/cst.cp*qt + dalphadT*cst.Ls/cst.cp*qt \
                + dalphadT*cst.Lv/cst.cp*qs - dalphadT*cst.Ls/cst.cp*qs \
                + alpha_w*cst.Lv/cst.cp*dqsatdT_w \
                + alpha_i*cst.Ls/cst.cp*dqsatdT_i

            tnr -= f / f_prime

        alpha_w = water_fraction(tnr)
        alpha_i = 1. - alpha_w
        qs = qsat(p, tnr)
        ql_qi = max(0., qt-qs)
        ql = alpha_w*ql_qi;
        qi = alpha_i*ql_qi;

        return tnr, ql, qi, qs