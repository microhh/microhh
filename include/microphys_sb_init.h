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

#include <iostream>

#include "constants.h"
#include "fast_math.h"
#include "microphys_sb06.h"

namespace Sb_init
{
    template<typename TF>
    inline TF gser(
            const TF a,
            const TF x)
    {
        const TF eps = 1e-7;
        const int itmax = 100;

        const TF gln = std::log(std::tgamma(a));

        if (x <= TF(0))
        {
            std::cout << "ERROR in GSER: x < 0" << std::endl;
            return TF(0);
        }

        TF ap = a;
        TF sum = TF(1) / a;
        TF del = sum;
        for (int n=0; n<itmax; ++n)
        {
            ap = ap + TF(1);
            del = del * x / ap;
            sum = sum + del;
            if (std::abs(del) < std::abs(sum) * eps)
                break;
        }

        if (std::abs(del) >= std::abs(sum)*eps)
        {
            std::cout << "ERROR in GSER: a too large, ITMAX too small" << std::endl;
            return TF(0);
        }

        return sum * std::exp(-x + a * std::log(x) - gln);
    };

    template<typename TF>
    inline TF gcf(
            const TF a,
            const TF x)
    {
        const int  itmax = 100;
        const TF eps = 3.e-7;
        const TF fpmin = 1.e-30;

        const TF gln = std::log(std::tgamma(a));

        TF b = x + TF(1) - a;
        TF c = TF(1) / fpmin;
        TF d = TF(1) / b;
        TF h = d;
        TF del;

        for (int i=0; i<itmax; ++i)
        {
            const TF an=-(i+1)*((i+1)-a);
            b = b + TF(2);
            d = an * d + b;

            if (std::abs(d) < fpmin)
                d = fpmin;

            c = b + an / c;

            if (std::abs(c) < fpmin)
                c = fpmin;

            d = TF(1) / d;
            del = d * c;
            h = h * del;

            if (std::abs(del-TF(1)) < eps)
                break;
        }

        if (std::abs(del-TF(1)) >= eps)
        {
            std::cout << "ERROR in GCF: a too large, ITMAX too small" << std::endl;
            return TF(0);
        }

        return std::exp(-x + a * std::log(x) - gln) * h;
    };

    template<typename TF>
    inline TF gammp(
            const TF a,
            const TF x)
    {
        if (x < TF(0) || a <= TF(0))
        {
            std::cout << "ERROR in GAMMP: bad arguments" << std::endl;
            return TF(0);
        }
        else if (x < a+1)
            return gser(a, x);
        else
        {
            const TF gammp = gcf(a, x);
            return TF(1) - gammp;
        }
    };

    template<typename TF>
    inline TF incgfct_lower(
            const TF a,
            const TF x)
    {
        const TF gam = gammp(a, x);
        const TF gln = std::log(std::tgamma(a));
        return std::exp(gln) * gam;
    };

    template<typename TF>
    void incgfct_lower_lookupcreate(
            const TF a,
            Gamlookuptable<TF>& ltable,
            const int nl,
            const int nlhr)
    {
        const TF c1 =  36.629433904824623;
        const TF c2 = -0.119475603955226;
        const TF c3 =  0.339332937820052;
        const TF c4 =  1.156369000458310;

        // Store parameters in the structure ltable:
        ltable.a = a;
        ltable.n = nl;

        // Allocate Memory for the table vectors:
        ltable.x.resize(nl);
        ltable.igf.resize(nl);

        // High-resolution part (not used, commented out for safety).
        //ltable.nhr = nlhr;
        //ltable.xhr.resize(nlhr);
        //ltable.igfhr.resize(nlhr);

        // Low resolution part of the table:
        // -----------------------------------
        // Maximum x-value of the lookup table (99.5-%-value):
        ltable.x[ltable.n-2] = c1 * ( TF(1) - std::exp(c2 * std::pow(a, c3)) ) + c4*a;

        // Create lookup table vectors:
        ltable.dx = ltable.x[ltable.n-2] / (ltable.n-TF(2));
        ltable.odx = TF(1) / ltable.dx;

        // Diese Schleife vektorisiert nicht wg. incgfct_lower():
        for (int i=0; i<ltable.n-1; i++)
        {
            ltable.x[i] = i * ltable.dx;
            ltable.igf[i] = incgfct_lower(a, ltable.x[i]);
        }

        // The last value is for x = infinity:
        ltable.x[ltable.n-1] = (ltable.n-1) * ltable.dx;
        ltable.igf[ltable.n-1] = std::tgamma(a);

        // High resolution part of the table (lowest 2% of the X-values):
        // -----------------------------------
        // NOTE BvS: for some reason, some of these functions (incgfct_lower...)
        // fail if you throw in a small dummy number for `nlhr`.
        // Create lookup table vectors:
        //const int index = std::round(0.01*(ltable.n-1)) - 1;
        //ltable.dxhr = ltable.x[index] / (ltable.nhr-TF(1));
        //ltable.odxhr = TF(1) / ltable.dxhr;

        //std::cout << index << " " << ltable.dxhr << std::endl;
        //throw 1;

        //// Diese Schleife vektorisiert nicht wg. incgfct_lower():
        //for (int i=0; i<ltable.nhr; ++i)
        //{
        //    ltable.xhr[i] = i * ltable.dxhr;
        //    ltable.igfhr[i] = incgfct_lower(a, ltable.xhr[i]);
        //}
    };



}
