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

#ifndef MICROPHYS_2MOM_WARM
#define MICROPHYS_2MOM_WARM

#include <cmath>

#include "microphys.h"

class Master;
class Input;
class Data_block;

template<typename> class Grid;
template<typename> class Stats;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Thermo;
template<typename> class Field3d;
template<typename> class Microphys;

// Constants specific or tuned for this microphysics scheme
namespace Micro_2mom_warm_constants
{
    template<typename TF> constexpr TF pi      = 3.14159265359;
    template<typename TF> constexpr TF Nc0     = 70e6;                // Fixed cloud droplet number
    template<typename TF> constexpr TF K_t     = 2.5e-2;              // Conductivity of heat [J/(sKm)]
    template<typename TF> constexpr TF D_v     = 3.e-5;               // Diffusivity of water vapor [m2/s]
    template<typename TF> constexpr TF rho_w   = 1.e3;                // Density water
    template<typename TF> constexpr TF rho_0   = 1.225;               // SB06, p48
    template<typename TF> constexpr TF pirhow  = pi<TF>*rho_w<TF>/6.;
    template<typename TF> constexpr TF mc_min  = 4.2e-15;             // Min mean mass of cloud droplet
    template<typename TF> constexpr TF mc_max  = 2.6e-10;             // Max mean mass of cloud droplet
    template<typename TF> constexpr TF mr_min  = mc_max<TF>;          // Min mean mass of precipitation drop
    template<typename TF> constexpr TF mr_max  = 3e-6;                // Max mean mass of precipitation drop // as in UCLA-LES
    template<typename TF> constexpr TF ql_min  = 1.e-6;               // Min cloud liquid water for which calculations are performed
    template<typename TF> constexpr TF qr_min  = 1.e-15;              // Min rain liquid water for which calculations are performed
}

template<typename TF>
class Microphys_2mom_warm : public Microphys<TF>
{
    public:
        Microphys_2mom_warm(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Microphys_2mom_warm();

        void init();
        void create(Input&, Data_block&, Stats<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, const double);

        void exec_stats(Stats<TF>&, std::string, Field3d<TF>&, Field3d<TF>&, Thermo<TF>&, const double);
        virtual void exec_dump(Dump<TF>&, unsigned long) {};
        virtual void exec_cross(Cross<TF>&, unsigned long);

        void get_mask(Field3d<TF>&, Field3d<TF>&, Stats<TF>&, std::string);
        bool has_mask(std::string);

        unsigned long get_time_limit(unsigned long, double);

        #ifdef USECUDA
        //void prepare_device() {};
        //void clear_device() {};
        //void forward_device() {};
        //void backward_device() {};
        #endif

    private:
        using Microphys<TF>::swmicrophys;
        using Microphys<TF>::master;
        using Microphys<TF>::grid;
        using Microphys<TF>::fields;

        Boundary_cyclic<TF> boundary_cyclic;

        bool swmicrobudget;     // Output full microphysics budget terms
        TF cflmax;              // Max CFL number in microphysics sedimentation

        std::vector<std::string> crosslist;                  // Cross-sections handled by this class
        std::vector<std::string> available_masks = {"qr"};   // Vector with the masks that fields can provide

        // Surface precipitation statistics
        std::vector<TF> rr_bot;   // 2D surface sedimentation flux (kg m-2 s-1 == mm s-1)
};
#endif
