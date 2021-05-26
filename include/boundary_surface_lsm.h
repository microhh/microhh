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

template<typename> class Diff;

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

    // Land surface
    std::vector<TF> rs;       // Surface resistance (canopy or soil, s m-1)
    std::vector<TF> H;        // Sensible heat flux (W m-2)
    std::vector<TF> LE;       // Latent heat flux (W m-2)
    std::vector<TF> G;        // Soil heat flux (W m-2)
    std::vector<TF> S;        // Storage flux (W m-2)
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
        void create(Input&, Netcdf_handle&, Stats<TF>&, Column<TF>&, Cross<TF>&);
        void set_values();

        const std::vector<TF>& get_z0m()  const { return z0m; };
        const std::vector<TF>& get_dudz() const { return dudz_mo; }
        const std::vector<TF>& get_dvdz() const { return dvdz_mo; }
        const std::vector<TF>& get_dbdz() const { return dbdz_mo; }

        void exec(Thermo<TF>&);
        void exec_stats(Stats<TF>&);
        void exec_column(Column<TF>&);
        void exec_cross(Cross<TF>&, unsigned long);

        void load(const int);
        void save(const int);

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();  // TMP BVS
        void backward_device(); // TMP BVS

        TF* get_z0m_g() { return z0m_g; };
        TF* get_ustar_g() { return ustar_g; };
        TF* get_obuk_g() { return obuk_g; };
        #endif

    protected:
        void process_input(Input&, Thermo<TF>&); // Process and check the surface input
        void init_surface(Input&); // Allocate and initialize the surface arrays
        void init_lsm(); // Allocate and initialize the surface arrays
        void init_solver(); // Prepare the lookup table's for the surface layer solver

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

        typedef std::map<std::string, Field3dBc<TF>> BcMap;
        using Boundary<TF>::sbc;

        bool sw_constant_z0;

        bool sw_homogeneous;
        bool sw_free_drainage;
        bool sw_water;
        bool sw_tile_stats;

        TF tskin_water;

        std::vector<std::string> tile_names {
            "veg", "soil" ,"wet"};
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
        int lookup_table_size;

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



        //#ifdef USECUDA
        //TF* z0m_g;
        //TF* z0h_g;
        //TF* obuk_g;
        //TF* ustar_g;
        //int* nobuk_g;

        //float* zL_sl_g;
        //float* f_sl_g;
        //#endif

        Boundary_type thermobc;

    protected:
        // Cross sections
        std::vector<std::string> cross_list;         // List of active cross variables

        void update_slave_bcs();
};
#endif
