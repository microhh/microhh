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

#ifndef DEPOSITION_H
#define DEPOSITION_H

#include <vector>
#include <string>
#include <map>
#include "timeloop.h"
#include "boundary_surface_lsm.h"
#include "boundary.h"


class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Boundary_surface_lsm;

/**
 * Class that creates a deposition liniked to the chemistry
 */

enum class Deposition_type {disabled, enabled, simple, average};

template<typename TF>
struct Deposition_tile
{
    std::string long_name;    // Descriptive name of tile

    std::vector<TF> vdo3;     // Deposition velocity of ozone (m s-1)
    std::vector<TF> vdno;     // Deposition velocity of no (m s-1)
    std::vector<TF> vdno2;    // Deposition velocity of no2 (m s-1)
    std::vector<TF> vdhno3;   // Deposition velocity of hno3 (m s-1)
    std::vector<TF> vdh2o2;   // Deposition velocity of h2o2 (m s-1)
    std::vector<TF> vdrooh;   // Deposition velocity of rooh (m s-1)
    std::vector<TF> vdhcho;   // Deposition velocity of hcho (m s-1)

    #ifdef USECUDA
    cuda_vector<TF> vdo3_g;
    cuda_vector<TF> vdno_g;
    cuda_vector<TF> vdno2_g;
    cuda_vector<TF> vdhno3_g;
    cuda_vector<TF> vdh2o2_g;
    cuda_vector<TF> vdrooh_g;
    cuda_vector<TF> vdhcho_g;
    #endif
};

template<typename TF>
using Deposition_tile_map = std::map<std::string, Deposition_tile<TF>>;

template<typename TF>
class Deposition
{
    public:
        Deposition(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the chemistry class.
        ~Deposition();                                       ///< Destructor  of the chemistry class.

        void init(Input&);                 ///< Initialize the arrays that contain the profiles.
        void create(Stats<TF>&, Cross<TF>&);
        void update_time_dependent(Timeloop<TF>&, Boundary<TF>&,
             TF*, TF*, TF*, TF*, TF*, TF*, TF*); ///< Update the time dependent deposition parameters.

        const TF get_vd(const std::string&) const;                  ///< get the standard vd value (o3, no, no2, ..)
        void get_tiled_mean(TF*, const std::string&, TF, const TF*, const TF*, const TF*);
        void update_vd_water(TF*, const std::string&, const TF*, const TF*, const int*, const TF*, const TF*);
        void exec_cross(Cross<TF>&, unsigned long);
        void spatial_avg_vd(TF*);

        #ifdef USECUDA
        void prepare_device();
        #endif

    protected:
        std::vector<std::string> cross_list;         // List of active cross variables

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_deposition;

        std::shared_ptr<Boundary_surface_lsm<TF>> boundary_surface_lsm;

        TF deposition_var;   // here put the local vars
        TF henry_so2;
        TF rsoil_so2;
        TF rwat_so2;
        TF rws_so2;

        std::vector<TF> rmes;
        std::vector<TF> rsoil;
        std::vector<TF> rcut;
        std::vector<TF> rws;
        std::vector<TF> rwat;
        std::vector<TF> diff;
        std::vector<TF> diff_scl;
        std::vector<TF> henry;
        std::vector<TF> f0;

        #ifdef USECUDA
        cuda_vector<TF> rmes_g;
        cuda_vector<TF> rsoil_g;
        cuda_vector<TF> rcut_g;
        cuda_vector<TF> rws_g;
        cuda_vector<TF> rwat_g;
        cuda_vector<TF> diff_g;
        cuda_vector<TF> diff_scl_g;
        cuda_vector<TF> henry_g;
        cuda_vector<TF> f0_g;
        #endif

        TF vd_o3;
        TF vd_no;
        TF vd_no2;
        TF vd_hno3;
        TF vd_h2o2;
        TF vd_rooh;
        TF vd_hcho;

        std::vector<std::string> deposition_tile_names {"veg", "soil" ,"wet"};
        Deposition_tile_map<TF> deposition_tiles;

        const std::string tend_name = "deposition";
        const std::string tend_longname = "Deposition";
};
#endif