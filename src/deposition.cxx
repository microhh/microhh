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

//#include <cstdio>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <utility>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"
#include "netcdf_interface.h"
#include "constants.h" 
#include "timeloop.h"
#include "deposition.h"
#include "boundary_surface_lsm.h"
#include "boundary.h"



namespace 
{


    template<typename TF>
    void calc_deposition(
	    const TF* vdo3,
	    const TF* vdno,
	    const TF* vdno2,
	    const TF* vdhno3,
	    const TF* vdh2o2,
	    const TF* vdrooh,
	    const TF* vdhcho,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk, const TF* const restrict dz, TF* const restrict rhoref)
    {
    TF o3 = 0.0;
    TF hcho = 0.0;
    int ndep = 0;
    for (int j=jstart; j<jend; ++j)
	for (int i=istart; i<iend; ++i)
	{
                const int ij = i + j*jj;
		o3   = o3   + vdo3[ij];
		hcho = hcho + vdhcho[ij];
		ndep = ndep + 1;
	}
    printf("mean o3 and hcho deposition %13.5e %13.5e m/s \n",o3/ndep,hcho/ndep);
    }
}

template<typename TF>
Deposition<TF>::Deposition(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
	
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();
    // Create surface tiles:
    for (auto& name : tile_names)
        tiles.emplace(name, Surface_tile<TF>{});
    tiles.at("veg" ).long_name = "vegetation";
    tiles.at("soil").long_name = "bare soil";
    tiles.at("wet" ).long_name = "wet skin";
}

template <typename TF>
Deposition<TF>::~Deposition()
{
}

template <typename TF>
void Deposition<TF>::init(Input& inputin)
{
    for (auto& it : fields.st)
    {
        const std::string type = inputin.get_item<std::string>("deposition", "swdeposition", it.first, "0");
        if (type == "0")
        {
            // Cycle to avoid reading unneeded namelist options.
            continue;
        }
        else if (type == "enabled")
        {
            cmap[it.first].type = Deposition_type::enabled;
        }
        else if (type == "simple")
        {
            cmap[it.first].type = Deposition_type::simple;
	}
        else if (type == "disabled")
        {
            cmap[it.first].type = Deposition_type::disabled;
        }
        else
            throw std::runtime_error("Invalid option for \"Deposition_type\"");
    }
    deposition_var = inputin.get_item<TF>("deposition", "deposition_var","", (TF)1e5);
    vd_o3   = inputin.get_item<TF>("deposition", "vdo3", "", (TF)0.005);
    vd_no   = inputin.get_item<TF>("deposition", "vdno", "", (TF)0.002);
    vd_no2  = inputin.get_item<TF>("deposition", "vdno2", "", (TF)0.005);
    vd_hno3 = inputin.get_item<TF>("deposition", "vdhno3", "", (TF)0.040);
    vd_h2o2 = inputin.get_item<TF>("deposition", "vdh2o2", "", (TF)0.018);
    vd_rooh = inputin.get_item<TF>("deposition", "vdrooh", "", (TF)0.008);
    vd_hcho = inputin.get_item<TF>("deposition", "vdhcho", "", (TF)0.0033);
}

template <typename TF>
void Deposition<TF>::create(const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    //

}


template <typename TF>
void Deposition<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary,
	    const TF* vdo3,
	    const TF* vdno,
	    const TF* vdno2,
	    const TF* vdhno3,
	    const TF* vdh2o2,
	    const TF* vdrooh,
	    const TF* vdhcho
		)
{
    for (auto& it : fields.st) if( cmap[it.first].type == Deposition_type::disabled) return;
    for (auto& it : fields.st) if( cmap[it.first].type == Deposition_type::simple) return;
    printf("Here is the entrypoint full deposition\n");
    auto& gd = grid.get_grid_data();
    const std::vector<TF>& z0m = boundary.get_z0m();
    printf("got z0m from boundary %13.3e  \n",z0m[20]);
    auto& tiles = boundary.get_tiles();

    calc_deposition<TF>(vdo3,vdno,vdno2,vdhno3,vdh2o2,vdrooh,vdhcho,
	    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
	    gd.icells, gd.ijcells, gd.dz.data(), fields.rhoref.data());

}

template<typename TF>
const TF Deposition<TF>::get_vd(const std::string& name) const
{
	if (name == "o3")
		return vd_o3;
	else if (name == "no")
		return vd_no;
	else if (name == "no2")
		return vd_no2;
	else if (name == "hno3")
		return vd_hno3;
	else if (name == "h2o2")
		return vd_h2o2;
	else if (name == "rooh")
		return vd_rooh;
	else if (name == "hcho")
		return vd_hcho;
        else  {
		std::string error = "Deposition::get_vd() can't return \"" + name + "\"";
		throw std::runtime_error(error);
	}
}

template class Deposition<double>;
//:template class Chemistry<float>;
