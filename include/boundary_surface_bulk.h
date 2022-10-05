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

#ifndef BOUNDARY_SURFACE_BULK_H
#define BOUNDARY_SURFACE_BULK_H

#include "boundary.h"
#include "stats.h"

template<typename> class Diff;

template<typename TF>
class Boundary_surface_bulk : public Boundary<TF>
{
    public:
        Boundary_surface_bulk(Master&, Grid<TF>&, Soil_grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_surface_bulk();

        void init(Input&, Thermo<TF>&);
        void create_cold_start(Netcdf_handle&);
        void create(Input&, Netcdf_handle&, Stats<TF>&, Column<TF>&, Cross<TF>&, Timeloop<TF>&);
        void set_values();

        const std::vector<TF>& get_z0m()  const { return z0m; };
        const std::vector<TF>& get_dudz() const { return dudz_mo; }
        const std::vector<TF>& get_dvdz() const { return dvdz_mo; }
        const std::vector<TF>& get_dbdz() const { return dbdz_mo; }

        void exec(Thermo<TF>&, Radiation<TF>&, Microphys<TF>&, Timeloop<TF>&);
        void exec_stats(Stats<TF>&);
        void exec_column(Column<TF>&);
        void exec_cross(Cross<TF>&, unsigned long) {};

        void load(const int, Thermo<TF>&);
        void save(const int, Thermo<TF>&);

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();  // TMP BVS
        void backward_device(); // TMP BVS

        TF* get_z0m_g()  { return z0m_g; };
        TF* get_dudz_g() { return dudz_mo_g; };
        TF* get_dvdz_g() { return dvdz_mo_g; };
        TF* get_dbdz_g() { return dbdz_mo_g; };
        #endif

    protected:
        void process_input(Input&, Thermo<TF>&); // Process and check the surface input
        void init_surface(Input&); // Allocate and initialize the surface arrays

    private:
        using Boundary<TF>::master;
        using Boundary<TF>::grid;
        using Boundary<TF>::soil_grid;
        using Boundary<TF>::fields;
        using Boundary<TF>::boundary_cyclic;
        using Boundary<TF>::swboundary;
        using Boundary<TF>::field3d_io;

        using Boundary<TF>::process_bcs;

        using Boundary<TF>::mbcbot;
        using Boundary<TF>::ubot;
        using Boundary<TF>::vbot;

        using Boundary<TF>::sbc;

        std::vector<TF> z0m;
        std::vector<TF> z0h;

        std::vector<TF> ustar;
        std::vector<TF> obuk;

        std::vector<TF> dudz_mo;
        std::vector<TF> dvdz_mo;
        std::vector<TF> dbdz_mo;

        #ifdef USECUDA
        TF* z0m_g;
        TF* obuk_g;

        TF* ustar_g;

        TF* dudz_mo_g;
        TF* dvdz_mo_g;
        TF* dbdz_mo_g;
        #endif

        // Transfer coefficients
        TF bulk_cm;
        std::map<std::string, TF> bulk_cs;




    protected:

        void update_slave_bcs();
};
#endif
