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
    void calc_tiled_mean(
            TF* const restrict fld,
            const TF* const restrict c_veg,
            const TF* const restrict c_soil,
            const TF* const restrict c_wet,
            const TF* const restrict fld_veg,
            const TF* const restrict fld_soil,
            const TF* const restrict fld_wet,
            const TF fac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                fld[ij] = (
                    c_veg [ij] * fld_veg [ij] +
                    c_soil[ij] * fld_soil[ij] +
                    c_wet [ij] * fld_wet [ij] ) * fac;
            }
    }
    template<typename TF>
    void calc_deposition_per_tile(
	    TF* restrict vdo3,
	    TF* restrict vdno,
	    TF* restrict vdno2,
	    TF* restrict vdhno3,
	    TF* restrict vdh2o2,
	    TF* restrict vdrooh,
	    TF* restrict vdhcho,
	    const TF* const restrict fraction,
            const int istart, const int iend, const int jstart, const int jend, const int jj)
    {
    TF o3 = 0.0;
    TF hcho = 0.0;
    TF frac = 0.0;
    int ndep = 0;
    for (int j=jstart; j<jend; ++j)
	for (int i=istart; i<iend; ++i)
	{
                const int ij = i + j*jj;
		ndep = ndep + 1;
		if (fraction[ij] < (TF) 1e-12) continue;
		frac = frac + fraction[ij];
		o3   = o3   + vdo3[ij]*fraction[ij];
		hcho = hcho + vdhcho[ij]*fraction[ij];
		
		//vdo3[ij] = 0.0;
		//vdno[ij] = 0.0;
		//vdno2[ij] = 0.0;
		//vdhno3[ij] = 0.0;
		//vdh2o2[ij] = 0.0;
		//vdrooh[ij] = 0.0;
		//vdhcho[ij] = 0.0;

	}
    if (frac > (TF) 1e-12) printf("mean fraction %13.5e vdo3 %13.5e vdhcho %13.5e \n",frac/ndep,o3/frac,hcho/frac);
    }
}

template<typename TF>
Deposition<TF>::Deposition(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
	
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();

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

    auto& gd = grid.get_grid_data();

    // Create surface tiles for deposition:
    for (auto& name : deposition_tile_names)
        deposition_tiles.emplace(name, Deposition_tile<TF>{});
    
    for (auto& tile : deposition_tiles){
        tile.second.vdo3.resize(gd.ijcells);
        tile.second.vdno.resize(gd.ijcells);
        tile.second.vdno2.resize(gd.ijcells);
        tile.second.vdno2.resize(gd.ijcells);
        tile.second.vdhno3.resize(gd.ijcells);
        tile.second.vdh2o2.resize(gd.ijcells);
        tile.second.vdrooh.resize(gd.ijcells);
        tile.second.vdhcho.resize(gd.ijcells);
    }

    deposition_tiles.at("veg" ).long_name = "vegetation";
    deposition_tiles.at("soil").long_name = "bare soil";
    deposition_tiles.at("wet" ).long_name = "wet skin";
    deposition_var = inputin.get_item<TF>("deposition", "deposition_var","", (TF)1e5);
    vd_o3   = inputin.get_item<TF>("deposition", "vdo3", "", (TF)0.005);
    vd_no   = inputin.get_item<TF>("deposition", "vdno", "", (TF)0.002);
    vd_no2  = inputin.get_item<TF>("deposition", "vdno2", "", (TF)0.005);
    vd_hno3 = inputin.get_item<TF>("deposition", "vdhno3", "", (TF)0.040);
    vd_h2o2 = inputin.get_item<TF>("deposition", "vdh2o2", "", (TF)0.018);
    vd_rooh = inputin.get_item<TF>("deposition", "vdrooh", "", (TF)0.008);
    vd_hcho = inputin.get_item<TF>("deposition", "vdhcho", "", (TF)0.0033);

    // fill also tile info
    for (auto& tile : deposition_tiles){
        std::fill(tile.second.vdo3.begin(),tile.second.vdo3.end(), vd_o3);
        std::fill(tile.second.vdno.begin(),tile.second.vdno.end(), vd_no);
        std::fill(tile.second.vdno2.begin(),tile.second.vdno2.end(), vd_no2);
        std::fill(tile.second.vdhno3.begin(),tile.second.vdhno3.end(), vd_hno3);
        std::fill(tile.second.vdh2o2.begin(),tile.second.vdh2o2.end(), vd_h2o2);
        std::fill(tile.second.vdrooh.begin(),tile.second.vdrooh.end(), vd_rooh);
        std::fill(tile.second.vdhcho.begin(),tile.second.vdhcho.end(), vd_hcho);
    }
}

template <typename TF>
void Deposition<TF>::create(const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc, Stats<TF>& stats)
{

}


template <typename TF>
void Deposition<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary,
	    TF* restrict vdo3,
	    TF* restrict vdno,
	    TF* restrict vdno2,
	    TF* restrict vdhno3,
	    TF* restrict vdh2o2,
	    TF* restrict vdrooh,
	    TF* restrict vdhcho
		)
{
    for (auto& it : fields.st) if( cmap[it.first].type == Deposition_type::disabled) return;
    for (auto& it : fields.st) if( cmap[it.first].type == Deposition_type::simple) return;
    auto& gd = grid.get_grid_data();

    // get information from lsm:
    auto& tiles = boundary.get_tiles();

    // calculate deposition per tile:
    for (auto& tile : tiles){
	    printf("Tiles from boundary_surface_lsm: %s  \n",tile.first.c_str());
	    calc_deposition_per_tile(
		    deposition_tiles.at(tile.first).vdo3.data(),
		    deposition_tiles.at(tile.first).vdno.data(),
		    deposition_tiles.at(tile.first).vdno2.data(),
		    deposition_tiles.at(tile.first).vdhno3.data(),
		    deposition_tiles.at(tile.first).vdh2o2.data(),
		    deposition_tiles.at(tile.first).vdrooh.data(),
		    deposition_tiles.at(tile.first).vdhcho.data(),
		    tile.second.fraction.data(),
	            gd.istart, gd.iend, gd.jstart, gd.jend, gd.icells);
    }
    // calculate tile-mean deposition for chemistry
    get_tiled_mean(vdo3,"o3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdno,"no",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdno2,"no2",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdhno3,"hno3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdh2o2,"h2o2",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdrooh,"rooh",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdhcho,"hcho",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
}

// routine to return standard deposition value:
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

template<typename TF>
void Deposition<TF>::get_tiled_mean(
    TF* restrict fld_out, std::string name, const TF fac, 
    const TF* const restrict fveg,
    const TF* const restrict fsoil,
    const TF* const restrict fwet)
{
    auto& gd = grid.get_grid_data();

    TF* fld_veg;
    TF* fld_soil;
    TF* fld_wet;

    // Yikes..
    if (name == "o3")
    {
        fld_veg  = deposition_tiles.at("veg").vdo3.data();
        fld_soil = deposition_tiles.at("soil").vdo3.data();
        fld_wet  = deposition_tiles.at("wet").vdo3.data();
    }
    else if (name == "no")
    {
        fld_veg  = deposition_tiles.at("veg").vdno.data();
        fld_soil = deposition_tiles.at("soil").vdno.data();
        fld_wet  = deposition_tiles.at("wet").vdno.data();
    }
    else if (name == "no2")
    {
        fld_veg  = deposition_tiles.at("veg").vdno2.data();
        fld_soil = deposition_tiles.at("soil").vdno2.data();
        fld_wet  = deposition_tiles.at("wet").vdno2.data();
    }
    else if (name == "hno3")
    {
        fld_veg  = deposition_tiles.at("veg").vdhno3.data();
        fld_soil = deposition_tiles.at("soil").vdhno3.data();
        fld_wet  = deposition_tiles.at("wet").vdhno3.data();
    }
    else if (name == "h2o2")
    {
        fld_veg  = deposition_tiles.at("veg").vdh2o2.data();
        fld_soil = deposition_tiles.at("soil").vdh2o2.data();
        fld_wet  = deposition_tiles.at("wet").vdh2o2.data();
    }
    else if (name == "rooh")
    {
        fld_veg  = deposition_tiles.at("veg").vdrooh.data();
        fld_soil = deposition_tiles.at("soil").vdrooh.data();
        fld_wet  = deposition_tiles.at("wet").vdrooh.data();
    }
    else if (name == "hcho")
    {
        fld_veg  = deposition_tiles.at("veg").vdhcho.data();
        fld_soil = deposition_tiles.at("soil").vdhcho.data();
        fld_wet  = deposition_tiles.at("wet").vdhcho.data();
    }
    else
        throw std::runtime_error("Cannot calculate tiled mean for variable \"" + name + "\"\\n");

    calc_tiled_mean(
            fld_out,
            fveg,
	    fsoil,
	    fwet,
            fld_veg,
            fld_soil,
            fld_wet,
            fac,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}

template class Deposition<double>;
//:template class Chemistry<float>;
