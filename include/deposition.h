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

/**
 * Class that creates a deposition liniked to the chemistry
 */

enum class Deposition_type {disabled, enabled, simple};

template<typename TF>
class Deposition
{
    public:
        Deposition(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the chemistry class.
        ~Deposition();                                       ///< Destructor  of the chemistry class.

        void init(Input&);                 ///< Initialize the arrays that contain the profiles.
        void create(const Timeloop<TF>&, std::string, Netcdf_handle&, Stats<TF>&);   ///< Read the profiles of the forces from the input.
        void update_time_dependent(Timeloop<TF>&, Boundary<TF>&,
			const TF*, const TF*, const TF*, const TF*, const TF*, const TF*, const TF*); ///< Update the time dependent deposition parameters.

        const TF get_vd(const std::string&) const;                  ///< get the standard vd value (o3, no, no2, ..)

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

	// internal variable
	struct Deposition_var
	{
	     Deposition_type type;
	};

        typedef std::map<std::string, Deposition_var> Deposition_map;

        Deposition_map cmap;

	TF deposition_var;   // here put the local vars

	TF vd_o3,vd_no,vd_no2,vd_hno3,vd_h2o2,vd_rooh,vd_hcho;

        std::vector<std::string> tile_names {"veg", "soil" ,"wet"};
        Tile_map<TF> tiles;


        const std::string tend_name = "deposition";
        const std::string tend_longname = "Deposition";
};
#endif
