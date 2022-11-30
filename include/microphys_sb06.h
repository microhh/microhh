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

#ifndef MICROPHYS_SB06_H
#define MICROPHYS_SB06_H

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

#include <cmath>
#include <map>

#include "microphys.h"
#include "field3d_operators.h"

class Master;
class Input;
class Netcdf_handle;

template<typename> class Grid;
template<typename> class Stats;
template<typename> class Dump;
template<typename> class Diff;
template<typename> class Cross;
template<typename> class Thermo;
template<typename> class Field3d;
template<typename> class Microphys;


// Derived type for hydrometeor species including pointers to data
template<typename TF>
struct Particle
{
    std::string name; // name of particle class
    TF nu;            // first shape parameter of size distribution
    TF mu;            // 2nd shape parameter
    TF x_max;         // max mean particle mass
    TF x_min;         // min mean particle mass
    TF a_geo;         // pre-factor in diameter-mass relation
    TF b_geo;         // exponent in diameter-mass relation
    TF a_vel;         // pre-factor in power law fall speed (all particles have a power law fall speed,
    TF b_vel;         // exponent in power law fall speed    some have also an Atlas-type relation)
    TF a_ven;         // first parameter in ventilation coefficient
    TF b_ven;         // 2nd parameter in ventilation coefficient
    TF cap;           // coefficient for capacity of particle
    TF vsedi_max;     // max bulk sedimentation velocity
    TF vsedi_min;     // min bulk sedimentation velocity
    TF* n;            // number density
    TF* q;            // mass density
    TF* rho_v;        // density correction of terminal fall velocity
};

// Because of OpenMP we have to separate the data pointers from the run-time-invariant coefficients.
// Therefore, we carry 2 data structures for each particle species, e.g. graupel and graupel_coeff.
// The following derived types are for the run-time coefficients
template<typename TF>
struct Particle_coeffs
{
    TF a_f; // ventilation coefficient, vent_coeff_a(particle,1)
    TF b_f; // ventilation coefficient, vent_coeff_b(particle,1) * N_sc**n_f / SQRT(nu_l)
    TF c_i; // 1.0/particle%cap
    TF c_z; // coefficient for 2nd mass moment
};

// .. non-spherical particles have an Atlas-type terminal fall velocity relation
template<typename TF>
struct Particle_nonsphere : public Particle_coeffs<TF>
{
    TF alfa; // 1st parameter in Atlas-type fall speed
    TF beta; // 2nd parameter in Atlas-type fall speed
    TF gama; // 3rd parameter in Atlas-type fall speed
};

// .. for spherical particles we need to store the coefficients for the
//    power law bulk sedimentation velocity
template<typename TF>
struct Particle_sphere : public Particle_coeffs<TF>
{
    TF coeff_alfa_n;
    TF coeff_alfa_q;
    TF coeff_lambda;
};

// raindrops have an Atlas-type terminal fall velocity relation
// and a mu-D-relation which is used in sedimentation and evaporation
// (see Seifert 2008, J. Atmos. Sci.)
template<typename TF>
struct Particle_rain_coeffs : public Particle_nonsphere<TF>
{
    TF cmu0; // Parameters for mu-D-relation of rain: max of left branch
    TF cmu1; // max of right branch
    TF cmu2; // scaling factor
    TF cmu3; // location of min value = breakup equilibrium diameter
    TF cmu4; // min value of relation
    TF cmu5; // exponent
};

template<typename TF>
struct Particle_cloud_coeffs : public Particle_nonsphere<TF>
{
    TF k_au; // ..Parameters for autoconversion
    TF k_sc; //     and selfcollection
};

template<typename TF>
struct Particle_ice_coeffs : public Particle_sphere<TF>
{
    TF sc_delta_n; //..Parameters for self-collection
    TF sc_delta_q; //   of cloud ice
    TF sc_theta_n;
    TF sc_theta_q;
};

template<typename TF>
struct Particle_snow_coeffs : public Particle_sphere<TF>
{
    TF sc_delta_n; //..Parameters for self-collection
    TF sc_theta_n; //   of snow
};

template<typename TF>
struct Particle_graupel_coeffs : public Particle_sphere<TF>
{
    TF sc_coll_n; //..Parameters for self-collection
};

template<typename TF>
struct Particle_frozen : public Particle<TF>
{
    TF ecoll_c;  //..maximum collision efficiency with cloud droplets
    TF D_crit_c; //..D-threshold for cloud riming
    TF q_crit_c; //..q-threshold for cloud riming
    TF s_vel;    //..dispersion of fall velocity for collection kernel (see SB2006, Eqs 60-63)
};


template<typename TF>
struct T_cfg_2mom
{
    int i2mom_solver;        // 0) explicit (1) semi-implicit solver
    int ccn_type;            // if not set by namelist, the ccn_type_gscp4 or ccn_type_gscp5 will win
    TF alpha_spacefilling;   // factor involved in the conversion of ice/snow to graupel by riming
    TF D_conv_ii;            // D-threshold for conversion to snow ice_selfcollection [m]
    TF D_rainfrz_ig;         // rain --> ice oder graupel [m]
    TF D_rainfrz_gh;         // rain --> graupel oder hail [m]
    TF rain_cmu0;            // asymptotic mue-value for small D_m in the mu-Dm-Relation of Seifert (2008)
    TF rain_cmu1;            // asymptotic mue-value for large D_m in the mu-Dm-Relation of Seifert (2008)
    TF rain_cmu3;            // D_br: equilibrium diameter for breakup and selfcollection, but only used in mue-D-relation [m]
    TF eva_q_fak_low;        // \  Parameters of the
    TF eva_q_fak_high;       // |  ramp-function eva_q_fak(D_m) for reduction
    TF eva_q_fak_Dbr_minfak; // |  of evaporation of drizzle-like rain
    TF eva_q_fak_Dbr_maxfak; // /
    TF melt_h_tune_fak;      // Factor to increase/decrease hail melting rate
    TF Tmax_gr_rime;         // Allow formation of graupel by riming ice/snow only at T < this threshold [K]
};

template<typename TF>
struct Collection_coeffs
{
    TF delta_n_aa;
    TF delta_n_ab;
    TF delta_n_bb;
    TF delta_q_aa;
    TF delta_q_ab;
    TF delta_q_bb;
    TF theta_n_aa;
    TF theta_n_ab;
    TF theta_n_bb;
    TF theta_q_aa;
    TF theta_q_ab;
    TF theta_q_bb;
};

//..these are coefficients for collection processes of the type a+b->c
template<typename TF>
struct Rain_riming_coeffs
{
    TF delta_n_aa;
    TF delta_n_ab;
    TF delta_n_bb;
    TF delta_q_aa;
    TF delta_q_ab;
    TF delta_q_ba;
    TF delta_q_bb;
    TF theta_n_aa;
    TF theta_n_ab;
    TF theta_n_bb;
    TF theta_q_aa;
    TF theta_q_ab;
    TF theta_q_ba;
    TF theta_q_bb;
};

template<typename TF>
struct Hydro_type
{
    std::string name;       // Species name (e.g. `rain`)
    std::string long_name;  // Long name (e.g. `rain specific humidity`)
    std::string units;      // Units (e.g. `kg kg-1`)
    bool is_mass;           // Switch between mass and density.

    std::vector<TF> precip_rate;

    // XY slices from tmp field, for implicit solver
    TF* v_sed_now;
    TF* v_sed_new;
    TF* flux_now;
    TF* flux_new;
    TF* sum;
    TF* impl;
    TF* slice;
    TF* conversion_tend;
    TF* tmp1;
};

template<typename TF>
class Microphys_sb06 : public Microphys<TF>
{
    public:
        Microphys_sb06(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Microphys_sb06();

        void init();
        void create(Input&, Netcdf_handle&, Stats<TF>&, Cross<TF>&, Dump<TF>&, Column<TF>&);
        void exec(Thermo<TF>&, const double, Stats<TF>&);

        void exec_stats(Stats<TF>&, Thermo<TF>&, const double);
        void exec_column(Column<TF>&);
        void exec_dump(Dump<TF>&, unsigned long) {};
        void exec_cross(Cross<TF>&, unsigned long);

        void get_mask(Stats<TF>&, std::string);
        bool has_mask(std::string);

        void get_surface_rain_rate(std::vector<TF>&);

        unsigned long get_time_limit(unsigned long, double);

        //#ifdef USECUDA
        //void get_surface_rain_rate_g(TF*);
        //void prepare_device();
        //void clear_device();
        //void backward_device();
        //#endif

    private:
        using Microphys<TF>::swmicrophys;
        using Microphys<TF>::master;
        using Microphys<TF>::grid;
        using Microphys<TF>::fields;
        using Microphys<TF>::field3d_operators;

        void init_2mom_scheme();
        void init_2mom_scheme_once();

        bool sw_warm;         // Switch between warm (true, old `2mom_warm`) and full (false) SB06 scheme
        bool sw_microbudget;  // Output full microphysics budget terms
        bool sw_debug;        // Print debug info.
        bool sw_integrate;    // Pre-integrate water species in implicit solver.
        double cfl_max;       // Max CFL number in microphysics sedimentation.

        const int cloud_type = 2673;
        const int mu_Dm_rain_typ = 1;
        TF rain_gfak;

        // Map with hydrometeor types.
        std::map<std::string, Hydro_type<TF>> hydro_types;

        // NOTE: this switch is set to True in ICON, but produces discontinuities in Nr and qr profiles.
        // Disable the feature for now, and discuss later with the ICON people.
        bool use_ql_sedi_rain = false;

        std::vector<std::string> crosslist; // Cross-sections handled by this class

        const std::string tend_name = "micro";
        const std::string tend_longname = "Microphysics";

        // Variables for microphysics.
        TF Nc0;  // Number concentration of cloud water (cm-3)
        TF Ni0;  // Number concentration of cloud ice (cm-3)

        //#ifdef USECUDA
        //TF* rr_bot_g;
        //TF* rs_bot_g;
        //TF* rg_bot_g;
        //#endif

        Particle<TF> cloud;
        Particle_cloud_coeffs<TF> cloud_coeffs;

        Particle<TF> rain;
        Particle_rain_coeffs<TF> rain_coeffs;

        Particle_frozen<TF> ice;
        Particle_ice_coeffs<TF> ice_coeffs;

        Particle_frozen<TF> snow;
        Particle_snow_coeffs<TF> snow_coeffs;

        Particle_frozen<TF> graupel;
        Particle_graupel_coeffs<TF> graupel_coeffs;

        Particle_frozen<TF> hail;
        Particle_sphere<TF> hail_coeffs;

        Collection_coeffs<TF> scr_coeffs;  // snow cloud riming
        Rain_riming_coeffs<TF> srr_coeffs; // snow rain riming
        Rain_riming_coeffs<TF> irr_coeffs; // ice rain riming
        Collection_coeffs<TF> icr_coeffs;  // ice cloud riming

        Collection_coeffs<TF> hrr_coeffs;  // hail rain riming
        Collection_coeffs<TF> grr_coeffs;  // graupel rain riming
        Collection_coeffs<TF> hcr_coeffs;  // hail cloud riming
        Collection_coeffs<TF> gcr_coeffs;  // graupel cloud  riming
        Collection_coeffs<TF> sic_coeffs;  // snow ice collection
        Collection_coeffs<TF> hic_coeffs;  // hail ice collection
        Collection_coeffs<TF> gic_coeffs;  // graupel ice collection
        Collection_coeffs<TF> hsc_coeffs;  // hail snow collection
        Collection_coeffs<TF> gsc_coeffs;  // graupel snow collection

        const Particle<TF> rainSBB = {
                "rainSBB", // name
                0.000000,  // nu
                0.333333,  // mu
                3.00E-06,  // x_max
                2.60E-10,  // x_min
                1.24E-01,  // a_geo
                0.333333,  // b_geo
                114.0137,  // a_vel
                0.234370,  // b_vel
                0.780000,  // a_ven
                0.308000,  // b_ven
                2.000000,  // cap
                2.000E+1,  // vsedi_max
                0.1,       // vsedi_min
                nullptr,   // n pointer
                nullptr,   // q pointer
                nullptr    // rho_v pointer
        };

        const Particle_rain_coeffs<TF> rainSBBcoeffs {
                0.0, 0.0, 0.0, 0.0,
                9.292000,  // ..alfa
                9.623000,  // ..beta
                6.222e+2,  // ..gama
                6.0000e0,  // ..cmu0
                3.000e+1,  // ..cmu1
                1.000e+3,  // ..cmu2
                1.100e-3,  // ..cmu3 = D_br
                1.0000e0,  // ..cmu4
                2          // ..cmu5
        };

        const Particle<TF> cloud_nue1mue1 {
                "cloud_nue1mue1", // name...Bezeichnung der Partikelklasse
                1.000000,         // nu.....Breiteparameter der Verteil.
                1.000000,         // mu.....Exp.-parameter der Verteil.
                2.60e-10,         // x_max..maximale Teilchenmasse D=80e-6m
                4.20e-15,         // x_min..minimale Teilchenmasse D=2.e-6m
                1.24e-01,         // a_geo..Koeff. Geometrie
                0.333333,         // b_geo..Koeff. Geometrie = 1/3
                3.75e+05,         // a_vel..Koeff. Fallgesetz
                0.666667,         // b_vel..Koeff. Fallgesetz
                0.780000,         // a_ven..Koeff. Ventilation (PK, S.541)
                0.308000,         // b_ven..Koeff. Ventilation (PK, S.541)
                2.00,             // cap....Koeff. Kapazitaet
                1.0,              // vsedi_max
                0.0,              // vsedi_min
                nullptr,          // n pointer
                nullptr,          // q pointer
                nullptr           // rho_v pointer
        };

        const Particle_frozen<TF> ice_cosmo5 {
                "ice_cosmo5", //.name...Bezeichnung der Partikelklasse
                0.000000,     //..nu...e..Breiteparameter der Verteil.
                0.333333,     //..mu.....Exp.-parameter der Verteil.
                1.00e-05,     //..x_max..maximale Teilchenmasse D=???e-2m
                1.00e-12,     //..x_min..minimale Teilchenmasse D=200e-6m
                0.835000,     //..a_geo..Koeff. Geometrie
                0.390000,     //..b_geo..Koeff. Geometrie
                2.77e+01,     //..a_vel..Koeff. Fallgesetz
                0.215790,     //..b_vel..Koeff. Fallgesetz = 0.41/1.9
                0.780000,     //..a_ven..Koeff. Ventilation (PK, S.541)
                0.308000,     //..b_ven..Koeff. Ventilation (PK, S.541)
                3.0,          //..cap....Koeff. Kapazitaet
                3.0,          //..vsedi_max
                0.0,          //..vsedi_min
                nullptr,      //..n pointer
                nullptr,      //..q pointer
                nullptr,      //..rho_v pointer
                0.80,         //..ecoll_c
                150.0e-6,     //..D_crit_c
                1.000e-5,     //..q_crit_c
                0.25          //..sigma_vel
        };

        const Particle_frozen<TF> snowSBB{
                "snowSBB", //..name...Bezeichnung der Partikelklasse
                0.000000,  //..nu.....Breiteparameter der Verteil.
                0.500000,  //..mu.....Exp.-parameter der Verteil.
                2.00e-05,  //..x_max..maximale Teilchenmasse
                1.00e-10,  //..x_min..minimale Teilchenmasse
                5.130000,  //..a_geo..Koeff. Geometrie, x = 0.038*D**2
                0.500000,  //..b_geo..Koeff. Geometrie = 1/2
                8.294000,  //..a_vel..Koeff. Fallgesetz
                0.125000,  //..b_vel..Koeff. Fallgesetz
                0.780000,  //..a_ven..Koeff. Ventilation (PK, S.541)
                0.308000,  //..b_ven..Koeff. Ventilation (PK, S.541)
                3.00,      //..cap....Koeff. Kapazitaet
                3.0,       //..vsedi_max
                0.1,       //..vsedi_min
                nullptr,   //..n pointer
                nullptr,   //..q pointer
                nullptr,   //..rho_v pointer
                0.80,      //..ecoll_c
                150.0e-6,  //..D_crit_c
                1.000e-5,  //..q_crit_c
                0.25       //..sigma_vel
        };

        const Particle_frozen<TF> graupelhail_cosmo5{
                "graupelhail_cosmo5", //.name
                1.000000,              //..nu.........1st shape parameter of the distribution
                0.333333,              //..mu.........2nd shape parameter of the distribution
                5.30e-04,              //..x_max......maximum particle mean mass
                4.19e-09,              //..x_min......minimum particle mean mass
                1.42e-01,              //..a_geo......particle geometry prefactor
                0.314000,              //..b_geo......particle geometry exponent = 1/3.10
                100.0,  //86.89371,    //..a_vel......terminal fall velocity prefactor
                0.34,   //0.268325,    //..b_vel......terminal fall velocity exponent
                0.780000,              //..a_ven......1st ventilation coefficient (PK, S.541)
                0.308000,              //..b_ven......2nd ventilation coefficient (PK, S.541)
                2.00,                  //..cap........capacity coefficient
                80.0,   //30.0,        //..vsedi_max..maximum bulk sedimentation velocity
                0.10,                  //..vsedi_min..minimum bulk sedimentation velocity
                nullptr,               //..n pointer..pointer to number density array
                nullptr,               //..q pointer..pointer to mass density array
                nullptr,               //..rho_v......pointer to density correction array
                1.0,                   //..ecoll_c....maximum collision efficiency with cloud droplets
                100.0e-6,              //..D_crit_c...D-threshold for cloud riming
                1.000e-6,              //..q_crit_c...q-threshold for cloud riming
                0.0                    //..sigma_vel..dispersion of fall velocity for collection kernel
        };

        const Particle_frozen<TF> hail_cosmo5{
                "hail_cosmo5", // !..name...Bezeichnung
                1.000000,      // !..nu.....Breiteparameter der Verteil.
                0.333333,      // !..mu.....Exp.-parameter der Verteil.
                5.00e-03,      // !..x_max..maximale Teilchenmasse
                2.60e-9,       // !..x_min..minimale Teilchenmasse
                0.1366,        // !..a_geo..Koeff. Geometrie
                0.333333,      // !..b_geo..Koeff. Geometrie = 1/3
                39.3,          // !..a_vel..Koeff. Fallgesetz
                0.166667,      // !..b_vel..Koeff. Fallgesetz
                0.780000,      // !..a_ven..Koeff. Ventilation (PK, S.541)
                0.308000,      // !..b_ven..Koeff. Ventilation (PK, S.541)
                2.00,          // !..cap....Koeff. Kapazitaet
                30.0,          // !..vsedi_max
                0.1,           // !..vsedi_min
                nullptr,       // !..n pointer
                nullptr,       // !..q pointer
                nullptr,       // !..rho_v pointer
                1.0,           // !..ecoll_c
                100.0e-6,      // !..D_crit_c
                1.000e-6,      // !..q_crit_c
                0.0            // !..sigma_vel
        };

        const T_cfg_2mom<TF> t_cfg_2mom{
                1,        // i2mom_solver: 0) explicit (1) semi-implicit solver
                -1,       // ccn_type: 6,7,8,9; if not set by namelist, the ccn_type_gscp4 or ccn_type_gscp5 will win
                0.01,     // alpha_spacefilling
                75.0e-6,  // D-threshold for conversion to snow ice_selfcollection [m]
                0.50e-3,  // D_rainfrz_ig [m]
                1.25e-3,  // D_rainfrz_gh [m]
                6.0,      // rain_cmu0
                30.0,     // rain_cmu1
                1.1e-3,   // rain_cmu3 = D_br (only in mue-D-relation, not in selfcollection-breakup-code!) [m]
                0.3,      // eva_q_fak_low
                1.0,      // eva_q_fak_high
                0.75,     // eva_q_fak_Dbr_minfak
                0.9,      // eva_q_fak_Dbr_maxfak
                1.0,      // melt_h_tune_fak
                270.16    // Tmax_gr_rime [K]
        };
};
#endif
