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

        bool swmicrobudget;     // Output full microphysics budget terms
        double cflmax;          // Max CFL number in microphysics sedimentation

        std::vector<std::string> crosslist; // Cross-sections handled by this class

        const std::string tend_name = "micro";
        const std::string tend_longname = "Microphysics";

        // Variables for microphysics.
        TF Nc0;  // Number concentration of cloud water (cm-3)

        std::vector<TF> rr_bot; // Rain rate at the bottom.
        std::vector<TF> rs_bot; // Snow rate at the bottom.
        std::vector<TF> rg_bot; // Graupel rate at the bottom.

        //#ifdef USECUDA
        //TF* rr_bot_g;
        //TF* rs_bot_g;
        //TF* rg_bot_g;
        //#endif

        Particle<TF> rain;
        Particle_rain_coeffs<TF> rain_coeffs;

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
};
#endif
