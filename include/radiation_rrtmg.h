/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 * Copyright (c) 2018-2019 Elynn Wu
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

#ifndef RADIATION_RRTMG
#define RADIATION_RRTMG

#include "radiation.h"
#include "field3d_operators.h"

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
class Radiation_rrtmg : public Radiation<TF>
{
	public:
		Radiation_rrtmg(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_rrtmg();

		bool check_field_exists(std::string name);
        void init();
        void create(Thermo<TF>&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&, Stats<TF>&);

		// Empty functions that should throw
		// Stats/dump not implemented in rrmtg now
        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&){throw 1;};

        void exec_stats(Stats<TF>&, Thermo<TF>&, Timeloop<TF>&){};
        void exec_cross(Cross<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&){};
        virtual void exec_dump(Dump<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&){};
        virtual void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&){};
	private:
		using Radiation<TF>::swradiation;
		using Radiation<TF>::master;
		using Radiation<TF>::grid;
		using Radiation<TF>::fields;
		using Radiation<TF>::field3d_operators;

		int ncol;
        int nlay;
        int nbndlw;

        int icld;
        int idrv;

        std::vector<double> play; // (ncol, nlay)
        std::vector<double> plev; // (ncol, nlay+1)
        std::vector<double> tlay; // (ncol, nlay)
        std::vector<double> tlev; // (ncol, nlay+1)

        std::vector<double> tsfc; // (ncol)

        std::vector<double> h2ovmr; // (ncol, nlay)
        std::vector<double> o3vmr;  // (ncol, nlay)
        std::vector<double> co2vmr; // (ncol, nlay)
        std::vector<double> ch4vmr; // (ncol, nlay)
        std::vector<double> n2ovmr; // (ncol, nlay)
        std::vector<double> o2vmr;  // (ncol, nlay)

        std::vector<double> cfc11vmr; // (ncol, nlay)
        std::vector<double> cfc12vmr; // (ncol, nlay)
        std::vector<double> cfc22vmr; // (ncol, nlay)
        std::vector<double> ccl4vmr;  // (ncol, nlay)
        std::vector<double> emis;     // (ncol, nbndlw)

        int inflglw;
        int iceflglw;
        int liqflglw;

        std::vector<double> cldfr;  // (ncol, nlay)
        std::vector<double> cicewp; // (ncol, nlay)
        std::vector<double> cliqwp; // (ncol, nlay)
        std::vector<double> reice;  // (ncol, nlay)
        std::vector<double> reliq;  // (ncol, nlay)
        std::vector<double> taucld; // (nbndlw, ncol, nlay)
        std::vector<double> tauaer; // (nbndlw, ncol, nlay)

        // OUTPUT
        std::vector<double> uflx;      // (ncol, nlay)
        std::vector<double> dflx;      // (ncol, nlay)
        std::vector<double> hr;        // (ncol, nlay)
        std::vector<double> uflxc;     // (ncol, nlay)
        std::vector<double> dflxc;     // (ncol, nlay)
        std::vector<double> hrc;       // (ncol, nlay)
        std::vector<double> duflx_dt;  // (ncol, nlay)
        std::vector<double> duflxc_dt; // (ncol, nlay)

        const std::string tend_name = "rad";
        const std::string tend_longname = "Radiation";
};
#endif
