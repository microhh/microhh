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

#ifndef MICROPHYS_NSW6_H
#define MICROPHYS_NSW6_H

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

template<typename TF>
class Microphys_nsw6 : public Microphys<TF>
{
    public:
        Microphys_nsw6(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Microphys_nsw6();

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
        double cflmax;          // Max CFL number in microphysics sedimentation

        std::vector<std::string> crosslist; // Cross-sections handled by this class

        const std::string tend_name = "micro";
        const std::string tend_longname = "Microphysics";

        // Variables for microphysics.
        TF Nc0; // Number concentration of cloud water (cm-3)

        std::vector<TF> rr_bot; // Rain rate at the bottom.
        std::vector<TF> rs_bot; // Snow rate at the bottom.
        std::vector<TF> rg_bot; // Graupel rate at the bottom.

        #ifdef USECUDA
        TF* rr_bot_g;
        TF* rs_bot_g;
        TF* rg_bot_g;
        #endif

};
#endif
