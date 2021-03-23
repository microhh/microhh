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

#ifndef BOUNDARY_SURFACE_TILED_H
#define BOUNDARY_SURFACE_TILED_H

#include "boundary.h"
#include "stats.h"

template<typename> class Diff;

template<typename TF>
struct MO_surface_tile
{
    std::vector<TF> obuk;   // Obukhov length (m)
    std::vector<TF> ustar;  // Friction velocity (m s-1)

    std::vector<TF> z0m;    // Roughness length momentum (m)
    std::vector<TF> z0h;    // Roughness length scalars (m)
};

template<typename TF>
using MO_tile_map = std::map<std::string, MO_surface_tile<TF>>;


template<typename TF>
class Boundary_surface_tiled : public Boundary<TF>
{
    public:
        Boundary_surface_tiled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_surface_tiled();

        void init(Input&, Thermo<TF>&);
        void create(Input&, Netcdf_handle&, Stats<TF>&, Column<TF>&, Cross<TF>&);
        void set_values();

        void get_ra(Field3d<TF>&);
        const std::vector<TF>& get_z0m() const { return z0m; };
        void get_duvdz(std::vector<TF>&, std::vector<TF>&);
        void get_dbdz(std::vector<TF>&, std::vector<TF>&);

        void calc_mo_stability(Thermo<TF>&, Land_surface<TF>&);
        void calc_mo_bcs_momentum(Thermo<TF>&, Land_surface<TF>&);
        void calc_mo_bcs_scalars(Thermo<TF>&, Land_surface<TF>&);

        void exec_stats(Stats<TF>&);
        void exec_column(Column<TF>&);
        void exec_cross(Cross<TF>&, unsigned long);

        void load(const int);
        void save(const int);

    protected:
        void check_settings(); // Check input settings
        void init_surface(Input&); // Allocate and initialize the surface arrays

    private:
        using Boundary<TF>::master;
        using Boundary<TF>::grid;
        using Boundary<TF>::fields;
        using Boundary<TF>::boundary_cyclic;
        using Boundary<TF>::swboundary;
        using Boundary<TF>::field3d_io;

        using Boundary<TF>::process_bcs;

        using Boundary<TF>::mbcbot;
        using Boundary<TF>::ubot;
        using Boundary<TF>::vbot;

        typedef std::map<std::string, Field3dBc<TF>> BcMap;
        using Boundary<TF>::sbc;

        // Switch between namelist or 2D field input z0m/h
        bool sw_constant_z0;

        // Tile properties
        MO_tile_map<TF> tiles;

        // Tile averaged quantities
        std::vector<TF> ustar;
        std::vector<TF> obuk;

        // Only used for `mlen` calculation in diffusion:
        std::vector<TF> z0m;

        Boundary_type thermobc;

    protected:
        // Cross sections
        std::vector<std::string> cross_list;         // List of active cross variables
        void update_slave_bcs();
};
#endif
