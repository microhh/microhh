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

#ifndef BOUNDARY_SURFACE_LSM_H
#define BOUNDARY_SURFACE_LSM_H

#include "boundary.h"
#include "stats.h"
#include "field3d_operators.h"

template<typename> class Diff;

enum class Soil_interpolation_type {Mean, Max, Harmonic_mean};

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

template<typename TF>
class Boundary_surface_lsm : public Boundary<TF>
{
    public:
        Boundary_surface_lsm(Master&, Grid<TF>&, Soil_grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_surface_lsm();

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
        void exec_cross(Cross<TF>&, unsigned long);

        void load(const int, Thermo<TF>&);
        void save(const int, Thermo<TF>&);

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();
        void backward_device();

        TF* get_z0m_g()  { return z0m_g; };
        TF* get_dudz_g() { return dudz_mo_g; };
        TF* get_dvdz_g() { return dvdz_mo_g; };
        TF* get_dbdz_g() { return dbdz_mo_g; };
        #endif

    protected:
        void process_input(Input&, Thermo<TF>&); // Process and check the surface input
        void init_surface_layer(Input&);         // Allocate and initialize the surface layer arrays
        void init_land_surface();                // Allocate and initialize the land surface arrays
        void init_solver();                      // Prepare the lookup table's for the surface layer solver
        void create_stats(Stats<TF>&, Column<TF>&, Cross<TF>&);

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

        Field3d_operators<TF> field3d_operators;

        typedef std::map<std::string, Field3dBc<TF>> BcMap;
        using Boundary<TF>::sbc;

        void get_tiled_mean(std::vector<TF>&, std::string, TF);

        bool sw_constant_z0;
        bool sw_homogeneous;
        bool sw_free_drainage;
        bool sw_water;
        bool sw_tile_stats;
        bool sw_tile_stats_col;
        bool sw_homogenize_sfc;

        TF emis_sfc;

        std::vector<std::string> tile_names {"veg", "soil" ,"wet"};
        Tile_map<TF> tiles;

        std::vector<float> zL_sl;
        std::vector<float> f_sl;
        std::vector<int> nobuk;

        std::vector<TF> z0m;
        std::vector<TF> z0h;

        std::vector<TF> ustar;
        std::vector<TF> obuk;

        std::vector<TF> dudz_mo;
        std::vector<TF> dvdz_mo;
        std::vector<TF> dbdz_mo;

        // Lookup tables van Genuchten parameterisation
        std::shared_ptr<Netcdf_file> nc_lookup_table;

        // Soil cross-sections
        std::vector<std::string> crosslist;

        // 2D fields with surface properties
        std::vector<TF> gD_coeff;        // Coefficient in response surface to VPD (Pa)
        std::vector<TF> c_veg;           // Vegetation fraction (-)
        std::vector<TF> lai;             // Leaf area index (-)
        std::vector<TF> rs_veg_min;      // Minimum vegetation resistance (s m-1)
        std::vector<TF> rs_soil_min;     // Minimum soil resistance (s m-1)
        std::vector<TF> lambda_stable;   // Skin conductivity stable conditions (W m-2 K-1)
        std::vector<TF> lambda_unstable; // Skin conductivity unstable conditions (W m-2 K-1)
        std::vector<TF> cs_veg;          // Heat capacity skin layer (J K-1 m-2)
        std::vector<int> water_mask;     // Mask for open water (-)
        std::vector<TF> t_bot_water;

        std::vector<TF> interception;   // Interception rain/dew by surface (m s-1)
        std::vector<TF> throughfall;    // Throughfall rain/dew onto soil (m s-1)
        std::vector<TF> infiltration;   // Infiltration moisture into soil (m s-1)
        std::vector<TF> runoff;         // Surface runoff from soil (m s-1)

        // Soil properties
        std::vector<int> soil_index;    // Index in lookup tables
        std::vector<TF> diffusivity;    // Full level (m2 s-1)
        std::vector<TF> diffusivity_h;  // Half level (m2 s-1)
        std::vector<TF> conductivity;   // Full level (unit m s-1)
        std::vector<TF> conductivity_h; // Half level (unit m s-1)
        std::vector<TF> source;         // Source term (unit s-1)
        std::vector<TF> root_fraction;  // Root fraction per soil layer (-)

        // Lookup table data obtained from input NetCDF file:
        std::vector<TF> theta_res;  // Residual soil moisture content (m3 m-3)
        std::vector<TF> theta_wp;   // Soil moisture content at wilting point (m3 m-3)
        std::vector<TF> theta_fc;   // Soil moisture content at field capacity (m3 m-3)
        std::vector<TF> theta_sat;  // Soil moisture content at saturation (m3 m-3)

        std::vector<TF> gamma_theta_sat;  // Conducticity soil moisture at saturation (m3 m-3)

        std::vector<TF> vg_a;  // van Genuchten parameter alpha (m-1)
        std::vector<TF> vg_l;  // van Genuchten parameter l (-)
        std::vector<TF> vg_n;  // van Genuchten parameter n (-)

        // Derived lookup table entries
        std::vector<TF> vg_m;  // van Genuchten parameter m (-)

        std::vector<TF> kappa_theta_max;  // Maximum diffusivity (m2 s-1)
        std::vector<TF> kappa_theta_min;  // Minimum diffusivity (m2 s-1)
        std::vector<TF> gamma_theta_max;  // Maximum conductivity (m s-1):
        std::vector<TF> gamma_theta_min;  // Minimum conductivity (m s-1)
        std::vector<TF> gamma_T_dry;      // Heat conductivity dry soil (m s-1)
        std::vector<TF> rho_C;            // Volumetric soil heat capacity (J m-3 K-1)

        #ifdef USECUDA
        void print_ij(const TF*);
        void get_tiled_mean_g(TF*, std::string, TF);

        // Surface layer:
        float* zL_sl_g;
        float* f_sl_g;
        int* nobuk_g;

        TF* z0m_g;
        TF* z0h_g;

        TF* ustar_g;
        TF* obuk_g;

        TF* dudz_mo_g;
        TF* dvdz_mo_g;
        TF* dbdz_mo_g;

        // Land-surface:
        TF* gD_coeff_g;        // Coefficient in response surface to VPD (Pa)
        TF* c_veg_g;           // Vegetation fraction (-)
        TF* lai_g;             // Leaf area index (-)
        TF* rs_veg_min_g;      // Minimum vegetation resistance (s m-1)
        TF* rs_soil_min_g;     // Minimum soil resistance (s m-1)
        TF* lambda_stable_g;   // Skin conductivity stable conditions (W m-2 K-1)
        TF* lambda_unstable_g; // Skin conductivity unstable conditions (W m-2 K-1)
        TF* cs_veg_g;          // Heat capacity skin layer (J K-1 m-2)
        int* water_mask_g;     // Mask for open water (-)
        TF* t_bot_water_g;

        TF* interception_g;   // Interception rain/dew by surface (m s-1)
        TF* throughfall_g;    // Throughfall rain/dew onto soil (m s-1)
        TF* infiltration_g;   // Infiltration moisture into soil (m s-1)
        TF* runoff_g;         // Surface runoff from soil (m s-1)

        // Soil properties
        int* soil_index_g;    // Index in lookup tables
        TF* diffusivity_g;    // Full level (m2 s-1)
        TF* diffusivity_h_g;  // Half level (m2 s-1)
        TF* conductivity_g;   // Full level (unit m s-1)
        TF* conductivity_h_g; // Half level (unit m s-1)
        TF* source_g;         // Source term (unit s-1)
        TF* root_fraction_g;  // Root fraction per soil layer (-)

        // Lookup table data obtained from input NetCDF file:
        TF* theta_res_g;  // Residual soil moisture content (m3 m-3)
        TF* theta_wp_g;   // Soil moisture content at wilting point (m3 m-3)
        TF* theta_fc_g;   // Soil moisture content at field capacity (m3 m-3)
        TF* theta_sat_g;  // Soil moisture content at saturation (m3 m-3)

        TF* gamma_theta_sat_g;  // Conducticity soil moisture at saturation (m3 m-3)

        TF* vg_a_g;  // van Genuchten parameter alpha (m-1)
        TF* vg_l_g;  // van Genuchten parameter l (-)
        TF* vg_n_g;  // van Genuchten parameter n (-)

        // Derived lookup table entries
        TF* vg_m_g;  // van Genuchten parameter m (-)

        TF* kappa_theta_max_g;  // Maximum diffusivity (m2 s-1)
        TF* kappa_theta_min_g;  // Minimum diffusivity (m2 s-1)
        TF* gamma_theta_max_g;  // Maximum conductivity (m s-1):
        TF* gamma_theta_min_g;  // Minimum conductivity (m s-1)
        TF* gamma_T_dry_g;      // Heat conductivity dry soil (m s-1)
        TF* rho_C_g;            // Volumetric soil heat capacity (J m-3 K-1)
        #endif

        Boundary_type thermobc;

    protected:
        // Cross sections
        std::vector<std::string> cross_list;         // List of active cross variables

        void update_slave_bcs();
};
#endif
