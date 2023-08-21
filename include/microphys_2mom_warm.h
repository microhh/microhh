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

#ifndef MICROPHYS_2MOM_WARM_H
#define MICROPHYS_2MOM_WARM_H

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

// Constants specific or tuned for this microphysics scheme
namespace Micro_2mom_warm_constants
{
    template<typename TF> constexpr TF pi       = 3.14159265359;
    template<typename TF> constexpr TF K_t      = 2.5e-2;              // Conductivity of heat [J/(sKm)]
    template<typename TF> constexpr TF D_v      = 3.e-5;               // Diffusivity of water vapor [m2/s]
    template<typename TF> constexpr TF rho_w    = 1.e3;                // Density water
    template<typename TF> constexpr TF rho_0    = 1.225;               // SB06, p48
    template<typename TF> constexpr TF pirhow   = pi<TF>*rho_w<TF>/6.;
    template<typename TF> constexpr TF mc_min   = 4.2e-15;             // Min mean mass of cloud droplet
    template<typename TF> constexpr TF mc_max   = 2.6e-10;             // Max mean mass of cloud droplet
    template<typename TF> constexpr TF mr_min   = mc_max<TF>;          // Min mean mass of precipitation drop
    template<typename TF> constexpr TF mr_max   = 3e-6;                // Max mean mass of precipitation drop // as in UCLA-LES
    template<typename TF> constexpr TF ql_min   = 1.e-6;               // Min cloud liquid water for which calculations are performed
    template<typename TF> constexpr TF qr_min   = 1.e-15;              // Min rain liquid water for which calculations are performed
    template<typename TF> constexpr TF cfl_min  = 1.e-5;               // Small non-zero limit at the CFL number
}

namespace Micro_2mom_warm_functions
{
    using namespace Micro_2mom_warm_constants;

    // Rational tanh approximation
    template<typename TF> CUDA_MACRO
    inline TF tanh2(const TF x)
    {
        return x * (TF(27) + x * x) / (TF(27) + TF(9) * x * x);
    }

    // Given rain water content (qr), number density (nr) and density (rho)
    // calculate mean mass of rain drop
    template<typename TF> CUDA_MACRO
    inline TF calc_rain_mass(const TF qr, const TF nr, const TF rho)
    {
        //TF mr = rho * qr / (nr + dsmall);
        #ifdef __CUDACC__
        TF mr = rho * qr / fmax(nr, TF(1.));
        return fmin(fmax(mr, mr_min<TF>), mr_max<TF>);
        #else
        TF mr = rho * qr / std::max(nr, TF(1.));
        return std::min(std::max(mr, mr_min<TF>), mr_max<TF>);
        #endif
    }

    // Given mean mass rain drop, calculate mean diameter
    template<typename TF> CUDA_MACRO
    inline TF calc_rain_diameter(const TF mr)
    {
        return pow(mr/pirhow<TF>, TF(1.)/TF(3.));
    }

    // Shape parameter mu_r
    template<typename TF> CUDA_MACRO
    inline TF calc_mu_r(const TF dr)
    {
        //return 1./3.; // SB06
        //return 10. * (1. + tanh(1200 * (dr - 0.0015))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
        return TF(10) * (TF(1) + tanh2(TF(1200) * (dr - TF(0.0015)))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
    }

    // Slope parameter lambda_r
    template<typename TF> CUDA_MACRO
    inline TF calc_lambda_r(const TF mur, const TF dr)
    {
        return pow((mur+3)*(mur+2)*(mur+1), TF(1.)/TF(3.)) / dr;
    }

    template<typename TF> CUDA_MACRO
    inline TF minmod(const TF a, const TF b)
    {
        #ifdef __CUDACC__
        return copysign(TF(1.), a) * fmax(TF(0.), fmin(fabs(a), TF(copysign(TF(1.), a))*b));
        #else
        return copysign(TF(1.), a) * std::max(TF(0.), std::min(std::abs(a), TF(copysign(TF(1.), a))*b));
        #endif
    }

    template<typename TF> CUDA_MACRO
    inline TF min3(const TF a, const TF b, const TF c)
    {
        #ifdef __CUDACC__
        return fmin(a, fmin(b, c));
        #else
        return std::min(a, std::min(b, c));
        #endif
    }

    template<typename TF> CUDA_MACRO
    inline TF max3(const TF a, const TF b, const TF c)
    {
        #ifdef __CUDACC__
        return fmax(a, fmax(b, c));
        #else
        return std::max(a, std::max(b, c));
        #endif
    }
}

template<typename TF>
class Microphys_2mom_warm : public Microphys<TF>
{
    public:
        Microphys_2mom_warm(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Microphys_2mom_warm();

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

        TF get_Nc0() { return this->Nc0; }
        TF get_Ni0() { return static_cast<TF>(1e5); } // CvH: this is a temporary fix with previous default value, Ni0 is 3D in tomita!

        unsigned long get_time_limit(unsigned long, double);

        #ifdef USECUDA
        void get_surface_rain_rate_g(TF*);
        void prepare_device();
        void clear_device();
        void backward_device();
        #endif

    private:
        using Microphys<TF>::swmicrophys;
        using Microphys<TF>::master;
        using Microphys<TF>::grid;
        using Microphys<TF>::fields;
        using Microphys<TF>::field3d_operators;

        bool swmicrobudget;     // Output full microphysics budget terms
        TF cflmax;              // Max CFL number in microphysics sedimentation

        std::vector<std::string> crosslist;                  // Cross-sections handled by this class
        std::vector<std::string> available_masks = {"qr"};   // Vector with the masks that fields can provide

        TF Nc0; // Cloud droplet number concentration.

        // Surface precipitation statistics
        std::vector<TF> rr_bot;   // 2D surface sedimentation flux (kg m-2 s-1 == mm s-1)

        #ifdef USECUDA
        TF* rr_bot_g;
        #endif

        const std::string tend_name = "micro";
        const std::string tend_longname = "Microphysics";
};
#endif
