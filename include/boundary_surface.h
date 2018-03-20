/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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

#ifndef BOUNDARY_SURFACE_H
#define BOUNDARY_SURFACE_H

#include "boundary.h"
#include "stats.h"

template<typename TF>
class Boundary_surface : public Boundary<TF>
{
    public:
        Boundary_surface(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_surface();

        void init(Input&, Thermo<TF>&);
        void create(Input&, Stats<TF>&);
        void set_values();

        void exec_stats(Mask<TF>&);
        void exec_cross(int);

        // Make these variables public for out-of-class usage.
        TF* obuk;
        int* nobuk;
        TF* ustar;

        TF z0m;
        TF z0h;

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();  // TMP BVS
        void backward_device(); // TMP BVS 

        TF* obuk_g;
        TF* ustar_g;
        int* nobuk_g;
        #endif

    protected:
        void process_input(Input&, Thermo<TF>&); // Process and check the surface input 
        void init_surface(); // Allocate and initialize the surface arrays
        void init_solver(); // Prepare the lookup table's for the surface layer solver
        void set_ustar(); // Set fixed ustar

    private:
        using Boundary<TF>::master;
        using Boundary<TF>::grid;
        using Boundary<TF>::fields;
        using Boundary<TF>::swboundary;

        using Boundary<TF>::process_bcs;

        using Boundary<TF>::mbcbot;
        using Boundary<TF>::ubot;
        using Boundary<TF>::vbot;

        typedef std::map<std::string, Field3dBc<TF>> BcMap;
        using Boundary<TF>::sbc;

        // surface scheme
        void update_bcs(Thermo<TF>&);

        void stability(TF*, TF*, TF*,
                       TF*, TF*, TF*,
                       TF*, TF*, TF*,
                       TF*, const TF*);
        void stability_neutral(TF*, TF*,
                               TF*, TF*,
                               TF*, TF*,
                               TF*, const TF*);
        void surfm(TF*, TF*,
                   TF*, TF*, TF*, TF*,
                   TF*, TF*, TF*, TF*,
                   TF, Boundary_type);
        void surfs(TF*, TF*, TF*,
                   TF*, TF*, TF*,
                   TF, Boundary_type);

        TF calc_obuk_noslip_flux     (const float* const, const float* const, int&, TF, TF, TF);
        TF calc_obuk_noslip_dirichlet(const float* const, const float* const, int&, TF, TF, TF);

        TF ustarin;

        float* zL_sl;
        float* f_sl;

        #ifdef USECUDA
        float* zL_sl_g;
        float* f_sl_g;
        #endif

        Boundary_type thermobc;

    protected:
        // cross sections
        // std::vector<std::string> crosslist;        // List with all crosses from ini file
        // std::vector<std::string> allowedcrossvars; // List with allowed cross variables

        void update_slave_bcs();
};
#endif
