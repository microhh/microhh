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

#ifndef RADIATION_PRESCRIBED_H
#define RADIATION_PRESCRIBED_H

#include "radiation.h"
#include "timedep.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Thermo;
template<typename> class Timeloop;
template<typename> class Timedep;

template<typename TF>
class Radiation_prescribed : public Radiation<TF>
{
	public:
		Radiation_prescribed(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_prescribed() {}

        void init(Timeloop<TF>&);
        void create(
                Input&, Netcdf_handle&, Thermo<TF>&,
                Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&, Stats<TF>&);

        unsigned long get_time_limit(unsigned long);
        std::vector<TF>& get_surface_radiation(std::string);
        void update_time_dependent(Timeloop<TF>&);

        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&)
        { throw std::runtime_error("\"get_radiation_field()\" is not implemented in radiation_rrtmpg"); }
		bool check_field_exists(std::string name)
        { throw std::runtime_error("\"check_field_exists()\" is not implemented in radiation_rrtmpg"); }

        // Empty functions which do nothing:
        void exec_all_stats(
                Stats<TF>&, Cross<TF>&, Dump<TF>&, Column<TF>&,
                Thermo<TF>&, Timeloop<TF>&, const unsigned long, const int) {};
        void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&) {};

	private:
		using Radiation<TF>::swradiation;
		using Radiation<TF>::master;
		using Radiation<TF>::grid;
		using Radiation<TF>::fields;
		using Radiation<TF>::field3d_operators;

        bool swtimedep_prescribed;

        // Input values
        TF lw_flux_dn;
        TF lw_flux_up;
        TF sw_flux_dn;
        TF sw_flux_up;

        // Surface radiative fluxes
        std::vector<TF> lw_flux_dn_sfc;
        std::vector<TF> lw_flux_up_sfc;

        std::vector<TF> sw_flux_dn_sfc;
        std::vector<TF> sw_flux_up_sfc;

        // Option for time varying surface fluxes
        std::unique_ptr<Timedep<TF>> tdep_sw_flux_dn;
        std::unique_ptr<Timedep<TF>> tdep_sw_flux_up;
        std::unique_ptr<Timedep<TF>> tdep_lw_flux_dn;
        std::unique_ptr<Timedep<TF>> tdep_lw_flux_up;
};
#endif
