/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#ifndef LAND_SURFACE
#define LAND_SURFACE

class Netcdf_file;

enum class Soil_interpolation_type {Mean, Max};

template<typename TF>
struct Surface_tile
{
    std::vector<TF> fraction;  // Grid point fraction tile (-)

    std::vector<TF> H;   // Sensible heat flux (W m-2)
    std::vector<TF> LE;  // Latent heat flux (W m-2)
    std::vector<TF> G;   // Soil heat flux (W m-2)

    std::vector<TF> T_bot;    // Skin temperature (K)
    std::vector<TF> qt_bot;   // Skin specific humidity (kg kg-1)
    std::vector<TF> thl_bot;  // Skin (liquid water) potential temperature (K)

    std::vector<TF> thl_fluxbot;  // Surface kinematic heat flux (K m s-1)
    std::vector<TF> qt_fluxbot;   // Surface kinematic moisture flux (kg kg-1 m s-1)
};

template<typename TF>
using Tile_map = std::map<std::string, Surface_tile<TF>>;

template<typename TF>
class Land_surface
{
    public:
        Land_surface(Master&, Grid<TF>&, Soil_grid<TF>&, Fields<TF>&, Input&);
        ~Land_surface();

        void init();
        void create_cold_start(Input&, Netcdf_handle&);
        void create_fields_grid_stats(Input&, Netcdf_handle&, Stats<TF>&, Cross<TF>&);

        void save_prognostic_fields(int);
        void load_prognostic_fields(int);

        void exec_soil();
        void exec_surface();
        void exec_stats(Stats<TF>&);
        void exec_cross(Cross<TF>&, unsigned long);

    private:
        Master& master;
        Grid<TF>& grid;
        Soil_grid<TF>& soil_grid;
        Fields<TF>& fields;

        bool sw_land_surface;
        bool sw_homogeneous;
        bool sw_free_drainage;

        // Land-surface properties
        Tile_map<TF> tiles;
        std::vector<TF> liquid_water_reservoir;  // Liquid water on leaves/surface (m)

        // Soil properties
        std::vector<int> soil_index;    // Index in lookup tables

        std::vector<TF> diffusivity;    // Full level (m2 s-1)
        std::vector<TF> diffusivity_h;  // Half level (m2 s-1)
        std::vector<TF> conductivity;   // Full level (unit m s-1)
        std::vector<TF> conductivity_h; // Half level (unit m s-1)
        std::vector<TF> source;         // Source term (unit s-1)

        // Soil cross-sections
        std::vector<std::string> crosslist;

        // Lookup tables van Genuchten parameterisation
        std::shared_ptr<Netcdf_file> nc_lookup_table;
        int lookup_table_size;

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
};
#endif
