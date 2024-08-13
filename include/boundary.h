/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <memory>

#include "timedep.h"
#include "boundary_cyclic.h"
#include "boundary_outflow.h"
#include "field3d_io.h"

class Master;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Soil_grid;
template<typename> class Fields;
template<typename> class Diff;
template<typename> class Thermo;
template<typename> class Timedep;
template<typename> class Stats;
template<typename> class Column;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Timeloop;
template<typename> class Radiation;
template<typename> class Microphys;

class Input;

enum class Boundary_type {Dirichlet_type, Neumann_type, Flux_type, Ustar_type, Off_type};
enum class Boundary_w_type {Normal_type, Conservation_type};
enum class Surface_model {Enabled, Disabled};

// Size of lookup table in Boundary_surface
const int nzL_lut = 10000;

/**
 * Structure containing the boundary options and values per 3d field.
 */
template<typename TF>
struct Field3dBc
{
    TF bot; ///< Value of the bottom boundary.
    TF top; ///< Value of the top boundary.
    Boundary_type bcbot; ///< Switch for the bottom boundary.
    Boundary_type bctop; ///< Switch for the top boundary.
};

template<typename TF>
struct Surface_tile
{
    std::string long_name;    // Descriptive name of tile

    // Shared
    std::vector<TF> fraction; // Grid point fraction tile (-)
    std::vector<TF> thl_bot;  // Skin (liquid water) potential temperature (K)
    std::vector<TF> qt_bot;   // Skin specific humidity (kg kg-1)

    // Surface layer
    std::vector<TF> obuk;     // Obukhov length (m)
    std::vector<TF> ustar;    // Friction velocity (m s-1)
    std::vector<TF> bfluxbot; // Friction velocity (m s-1)
    std::vector<int> nobuk;   // Index in LUT
    std::vector<TF> ra;       // Aerodynamic resistance (s m-1)

    // Land surface
    std::vector<TF> rs;       // Surface resistance (canopy or soil, s m-1)
    std::vector<TF> H;        // Sensible heat flux (W m-2)
    std::vector<TF> LE;       // Latent heat flux (W m-2)
    std::vector<TF> G;        // Soil heat flux (W m-2)
    std::vector<TF> S;        // Storage flux (W m-2)

    #ifdef USECUDA
    // Shared
    TF* fraction_g; // Grid point fraction tile (-)
    TF* thl_bot_g;  // Skin (liquid water) potential temperature (K)
    TF* qt_bot_g;   // Skin specific humidity (kg kg-1)

    // Surface layer
    TF* obuk_g;     // Obukhov length (m)
    TF* ustar_g;    // Friction velocity (m s-1)
    TF* bfluxbot_g; // Friction velocity (m s-1)
    int* nobuk_g;   // Index in LUT
    TF* ra_g;       // Aerodynamic resistance (s m-1)

    // Land surface
    TF* rs_g;       // Surface resistance (canopy or soil, s m-1)
    TF* H_g;        // Sensible heat flux (W m-2)
    TF* LE_g;       // Latent heat flux (W m-2)
    TF* G_g;        // Soil heat flux (W m-2)
    TF* S_g;        // Storage flux (W m-2)
    #endif
};

template<typename TF>
using Tile_map = std::map<std::string, Surface_tile<TF>>;

/**
 * Base class for the boundary scheme.
 * This class handles the case when the boundary is turned off. Derived classes are
 * implemented that handle different boundary schemes.
 */
template<typename TF>
class Boundary
{
    public:
        Boundary(Master&, Grid<TF>&, Soil_grid<TF>&, Fields<TF>&, Input&); ///< Constuctor of the boundary class.
        virtual ~Boundary(); ///< Destructor of the boundary class.

        static std::shared_ptr<Boundary> factory(
            Master&, Grid<TF>&, Soil_grid<TF>&, Fields<TF>&, Input&); ///< Factory function for boundary class generation.

        virtual void init(Input&, Thermo<TF>&, const Sim_mode);   ///< Initialize the fields.
        virtual void create_cold_start(Netcdf_handle&); ///< Create fields for cold start.
        virtual void create(
                Input&, Netcdf_handle&, Stats<TF>&, Column<TF>&,
                Cross<TF>&, Timeloop<TF>&); ///< Create the fields.

        virtual void update_time_dependent(Timeloop<TF>&); ///< Update the time dependent parameters.

        virtual void set_values(); ///< Set all 2d fields to the prober BC value.

        virtual void set_ghost_cells(); ///< Set the top and bottom ghost cells
        virtual void set_ghost_cells_w(Boundary_w_type); ///< Update the boundary conditions.

        void set_prognostic_cyclic_bcs();
        void set_prognostic_outflow_bcs();

        virtual void exec(Thermo<TF>&, Radiation<TF>&, Microphys<TF>&, Timeloop<TF>&);
        virtual void exec_stats(Stats<TF>&); ///< Execute statistics of surface
        virtual void exec_column(Column<TF>&); ///< Execute column statistics of surface
        virtual void exec_cross(Cross<TF>&, unsigned long) {}; ///< Execute cross statistics of surface

        virtual void load(const int, Thermo<TF>&) {};
        virtual void save(const int, Thermo<TF>&) {};

        // Get functions for various 2D fields
        virtual const std::vector<TF>& get_z0m() const;
        virtual const std::vector<TF>& get_dudz() const;
        virtual const std::vector<TF>& get_dvdz() const;
        virtual const std::vector<TF>& get_dbdz() const;
        virtual const std::vector<TF>& get_lai() const;
        virtual const std::vector<int>& get_water_mask() const;
        virtual const std::vector<TF>& get_c_veg() const;
        virtual const Tile_map<TF>& get_tiles() const;

        std::string get_switch();

        #ifdef USECUDA
        virtual cuda_vector<TF>& get_z0m_g();
        virtual cuda_vector<TF>& get_dudz_g();
        virtual cuda_vector<TF>& get_dvdz_g();
        virtual cuda_vector<TF>& get_dbdz_g();

        virtual void prepare_device(Thermo<TF>&);
        virtual void forward_device(Thermo<TF>&);
        virtual void backward_device(Thermo<TF>&);
        virtual void clear_device(Thermo<TF>&);
        #endif

    protected:
        Master& master;
        Grid<TF>& grid;
        Soil_grid<TF>& soil_grid;
        Fields<TF>& fields;
        Boundary_cyclic<TF> boundary_cyclic;
        Field3d_io<TF> field3d_io;
        Boundary_outflow<TF> boundary_outflow;

        std::string swboundary;

        Boundary_type mbcbot;
        Boundary_type mbctop;

        TF ubot;
        TF utop;
        TF vbot;
        TF vtop;

        typedef std::map<std::string, Field3dBc<TF>> BcMap;
        BcMap sbc;

        std::map<std::string, Timedep<TF>*> tdep_bc;

        // Spatial sbot input:
        std::vector<std::string> sbot_2d_list;

        // Scalar in/outflow
        std::vector<std::string> scalar_outflow;
        std::map<std::string, std::vector<TF>> inflow_profiles;
        bool swtimedep_outflow;
        std::map<std::string, Timedep<TF>*> tdep_outflow;

        // Time varying spatial sbot input:

        bool swtimedep_sbot_2d;
        std::map<std::string, std::vector<TF>> sbot_2d_prev;
        std::map<std::string, std::vector<TF>> sbot_2d_next;
        unsigned long iloadtime_sbot_2d;
        unsigned long itime_sbot_2d_prev;
        unsigned long itime_sbot_2d_next;

        void process_bcs(Input&); ///< Process the boundary condition settings from the ini file.
        void process_time_dependent(Input&, Netcdf_handle&, Timeloop<TF>&); ///< Process the time dependent settings from the ini file.
        void process_inflow(Input&, Netcdf_handle&); ///< Process the time dependent settings from the ini file.

        #ifdef USECUDA
        std::map<std::string, TF*> inflow_profiles_g;
        std::map<std::string, TF*> sbot_2d_prev_g;
        std::map<std::string, TF*> sbot_2d_next_g;
        #endif

        // GPU functions and variables
        void set_bc_g   (TF*, TF*, TF*, Boundary_type, TF, TF, TF, bool); ///< Set the values for the boundary fields.
        void set_bc_2d_g(TF*, TF*, TF*, TF*, Boundary_type, TF, TF, bool); ///< Set the values for the boundary fields.

    private:
        virtual void update_slave_bcs(); // Update the slave boundary values.
};
#endif
