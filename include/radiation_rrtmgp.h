/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#ifndef RADIATION_RRTMGP_H
#define RADIATION_RRTMGP_H

#include "radiation.h"
#include "field3d_operators.h"

#include "Gas_concs.h"
#include "Gas_optics.h"
#include "Source_functions.h"
#include "Cloud_optics.h"

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

template<typename TF>
class Radiation_rrtmgp : public Radiation<TF>
{
	public:
		Radiation_rrtmgp(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_rrtmgp() {}

		bool check_field_exists(std::string name)
        { throw std::runtime_error("Not implemented"); }

        void init();
        void create(
                Input&, Netcdf_handle&, Thermo<TF>&,
                Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&, Stats<TF>&);

        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&)
        { throw std::runtime_error("Not implemented"); }

        void exec_all_stats(
                Stats<TF>&, Cross<TF>&, Dump<TF>&,
                Thermo<TF>&, Timeloop<TF>&,
                const unsigned long, const int);

        void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&) {};

	private:
		using Radiation<TF>::swradiation;
		using Radiation<TF>::master;
		using Radiation<TF>::grid;
		using Radiation<TF>::fields;
		using Radiation<TF>::field3d_operators;

        void create_column(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&);
        void create_column_longwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&,
                const Gas_concs<double>&);
        void create_column_shortwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&,
                const Gas_concs<double>&);

        void create_solver(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&);
        void create_solver_longwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&,
                const Gas_concs<double>&);
        void create_solver_shortwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&,
                const Gas_concs<double>&);

        void exec_longwave(
                Thermo<TF>&, Timeloop<TF>&, Stats<TF>&,
                Array<double,2>&, Array<double,2>&, Array<double,2>&,
                const Array<double,2>&, const Array<double,2>&,
                const Array<double,2>&, const Array<double,2>&);

        void exec_shortwave(
                Thermo<TF>&, Timeloop<TF>&, Stats<TF>&,
                Array<double,2>&, Array<double,2>&, Array<double,2>&, Array<double,2>&,
                const Array<double,2>&, const Array<double,2>&,
                const Array<double,2>&, const Array<double,2>&);

        void exec_stats(Stats<TF>&, Thermo<TF>&, Timeloop<TF>&);
        void exec_cross(Cross<TF>&, const int, Thermo<TF>&, Timeloop<TF>&);
        void exec_dump(Dump<TF>&, const int, Thermo<TF>&, Timeloop<TF>&) {};

        const std::string tend_name = "rad";
        const std::string tend_longname = "Radiation";

        bool sw_longwave;
        bool sw_shortwave;
        double dt_rad;
        double next_rad_time;

        std::vector<std::string> crosslist;

        // RRTMGP related variables.
        double tsi_scaling;
        double t_sfc;
        double emis_sfc;
        double sfc_alb_dir;
        double sfc_alb_dif;
        double mu0;

        // The reference column for the full profile.
        Array<double,2> lw_flux_dn_inc;
        Array<double,2> sw_flux_dn_dir_inc;
        Array<double,2> sw_flux_dn_dif_inc;

        // The full solver.
        Gas_concs<double> gas_concs;
        std::unique_ptr<Gas_optics<double>> kdist_lw;
        std::unique_ptr<Gas_optics<double>> kdist_sw;

        std::unique_ptr<Cloud_optics<double>> cloud_lw;
        std::unique_ptr<Cloud_optics<double>> cloud_sw;
};
#endif
