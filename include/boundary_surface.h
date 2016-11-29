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

#ifndef BOUNDARY_SURFACE
#define BOUNDARY_SURFACE

#include "boundary.h"
#include "stats.h"

class Model;
class Input;
class Stats;
struct Mask;

class Boundary_surface : public Boundary
{
    public:
        Boundary_surface(Model*, Input*);
        ~Boundary_surface();

        virtual void init(Input*);
        void create(Input*);
        virtual void set_values();

        void exec_stats(Mask*); ///< Execute statistics of surface
        void exec_cross();      ///< Execute cross sections of surface

        // Make these variables public for out-of-class usage.
        double* obuk;
        int*    nobuk;
        double* ustar;

        double z0m;
        double z0h;

#ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();  // TMP BVS
        void backward_device(); // TMP BVS 

        double* obuk_g;
        double* ustar_g;
        int*    nobuk_g;
#endif

    protected:
        void process_input(Input *);   // Process and check the surface input 
        void init_surface();           // Allocate and initialize the surface arrays
        void init_solver();            // Prepare the lookup table's for the surface layer solver
        void set_ustar();              // Set fixed ustar

    private:

        // surface scheme
        void update_bcs();

        void stability(double*, double*, double*,
                       double*, double*, double*,
                       double*, double*, double*,
                       double*, double*);
        void stability_neutral(double*, double*,
                               double*, double*,
                               double*, double*,
                               double*, double*);
        void surfm(double*, double*,
                   double*, double*, double*, double*,
                   double*, double*, double*, double*,
                   double, int);
        void surfs(double*, double*, double*,
                   double*, double*, double*,
                   double, int);

        double calc_obuk_noslip_flux     (const float* const, const float* const, int&, double, double, double);
        double calc_obuk_noslip_dirichlet(const float* const, const float* const, int&, double, double, double);

        double ustarin;


        float* zL_sl;
        float* f_sl;

#ifdef USECUDA
        float* zL_sl_g;
        float* f_sl_g;
#endif
        int thermobc;

    protected:
        // cross sections
        std::vector<std::string> crosslist;        // List with all crosses from ini file
        std::vector<std::string> allowedcrossvars; // List with allowed cross variables

        Stats* stats;
        void update_slave_bcs();
};
#endif
