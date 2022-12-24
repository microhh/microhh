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
#include "netcdf_interface.h"

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
            //std::cout << "ERROR in GSER: x < 0" << std::endl;
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
            //std::cout << "ERROR in GSER: a too large, ITMAX too small" << std::endl;
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


    template<typename TF>
    void init_dmin_wg_gr_ltab_equi(
            Lookupt_4D<TF>& ltab,
            const std::string& file_name,
            Particle<TF>& particle,
            Master& master)
    {
        /*
            Function for setting up the wet growth diameter of a frozen hydrometeor
            type as function of supercooled LWC qw, frozen content qi, pressure p and temperature T.
            Is needed for the Parameterization of conversion from graupel to hail
            via wet growth of graupel.

            A corresponding 4D lookup table is read from an external file and is made
            equidistant along all table dimensions for better vectorization of table lookup
            (quadro-linear interpolation). Here, 3 of the dimensions (qw, qi, p) are already assumed
            to be equidistant in the table file. Only T can be non-equidistant and is
            made equidistant by linear interpolation.

            NOTE BvS: I interpolated the table offline (which is a few lines of Python code with Xarray);
            the interpolation functionality from ICON is therefore not available in this function!!
        */

        master.print_message("Reading lookup table \"%s\" for %s\n", file_name.c_str(), particle.name.c_str());

        // Read the lookup table.
        Netcdf_file coef_nc(master, file_name, Netcdf_mode::Read);

        // NOTE to self: indexing is [n4,n3,n2,n1] = [qi,qw,T,p]

        // Set dimensions in lookup table.
        ltab.n4 = coef_nc.get_dimension_size("nqi");
        ltab.n3 = coef_nc.get_dimension_size("nqw");
        ltab.n2 = coef_nc.get_dimension_size("ntemp");
        ltab.n1 = coef_nc.get_dimension_size("npres");
        ltab.n  = ltab.n1 * ltab.n2 * ltab.n3 * ltab.n4;

        // Set strides in lookup table for 1D access in 4D array.
        ltab.n1_stride = 1;
        ltab.n2_stride = ltab.n1;
        ltab.n3_stride = ltab.n1 * ltab.n2;
        ltab.n4_stride = ltab.n1 * ltab.n2 * ltab.n3;

        // Check global attributes (a_geo et al.) with values in `particle`.
        bool success = true;
        auto check = [&](const std::string& name, const TF ref_value)
        {
            const TF nc_value = coef_nc.get_attribute<TF>(name, "global");
            if (std::abs(ref_value - nc_value) > 1e-12)
            {
                master.print_warning(
                    "Wrong \"%s\" in table; value from table=%f, from particle=%f\n", name.c_str(), nc_value, ref_value);
                success=false;
            }
        };

        check("a_geo", particle.a_geo);
        check("b_geo", particle.b_geo);
        check("a_vel", particle.a_vel);
        check("b_vel", particle.b_vel);

        if (success)
            master.print_message("Check: a_geo/b_geo/a_vel/b_vel values are correct!\n");
        else
            throw std::runtime_error("One of more coefficients incorrect in lookup table.");

        // Allocate vectors.
        ltab.x1.resize(ltab.n1);
        ltab.x2.resize(ltab.n2);
        ltab.x3.resize(ltab.n3);
        ltab.x4.resize(ltab.n4);
        ltab.ltable.resize(ltab.n);

        // Read data from NetCDF.
        coef_nc.get_variable(ltab.x1, "p", {0}, {ltab.n1});
        coef_nc.get_variable(ltab.x2, "T", {0}, {ltab.n2});
        coef_nc.get_variable(ltab.x3, "qw", {0}, {ltab.n3});
        coef_nc.get_variable(ltab.x4, "qi", {0}, {ltab.n4});
        coef_nc.get_variable(ltab.ltable, "Dmin_wetgrowth_table",
                {0, 0, 0, 0}, {ltab.n4, ltab.n3, ltab.n2, ltab.n1});

        // Get fill_value from Dmin_wetgrowth_table, and mask lookup table.
        const TF dmin_fillval = coef_nc.get_attribute<TF>("_FillValue", "Dmin_wetgrowth_table");

        for (int i=0; i<ltab.n; ++i)
        {
            if (std::abs(ltab.ltable[i] - dmin_fillval) < 1e-6)
                ltab.ltable[i] = TF(999.99);
        }

        // Unit conversion NetCDF file.
        for (int i=0; i<ltab.n1; ++i)
            ltab.x1[i] *= TF(100);                  // Conversion from hPa to P

        for (int i=0; i<ltab.n2; ++i)
            ltab.x2[i] += Constants::T0<TF>;        // Conversion from deg C to K

        for (int i=0; i<ltab.n3; ++i)
            ltab.x3[i] *= TF(0.001);                // Conversion from g/m3 to kg/m3

        for (int i=0; i<ltab.n4; ++i)
            ltab.x4[i] *= TF(0.001);                // Conversion from g/m3 to kg/m3

        for (int i=0; i<ltab.n; ++i)
            ltab.ltable[i] *= TF(0.001);           // Conversion from mm to m.

        // Set dimension spacings.
        ltab.dx1   = ltab.x1[1] - ltab.x1[0];
        ltab.odx1  = TF(1) / ltab.dx1;

        ltab.dx2   = ltab.x2[1] - ltab.x2[0];
        ltab.odx2  = TF(1) / ltab.dx2;

        ltab.dx3   = ltab.x3[1] - ltab.x3[0];
        ltab.odx3  = TF(1) / ltab.dx3;

        ltab.dx4   = ltab.x4[1] - ltab.x4[0];
        ltab.odx4  = TF(1) / ltab.dx4;

        // Check to make sure that users are not accidently using the non-interpolated version.
        for (int i=0; i<ltab.n2-1; ++i)
        {
            if (std::abs((ltab.x2[i+1] - ltab.x2[i]) - ltab.dx2) > 1e-12)
                throw std::runtime_error("Lookup table not equidistant in T-dimension!");
        }

        master.print_message("Table \"%s\" parsed successfully!\n", file_name.c_str());




    };
}
