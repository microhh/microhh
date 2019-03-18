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
#ifndef RADIATION_GCSS
#define RADIATION_GCSS
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <math.h> /* floor */
#include <time.h>
#include <vector>

#include "field3d_operators.h"
#include "radiation.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "data_block.h"
#include "stats.h"
#include "cross.h"
#include "dump.h"
#include "column.h"
#include "constants.h"
#include "timeloop.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Timedep;
template<typename> class Timeloop;

template<typename TF>
class Radiation_gcss : public Radiation<TF>
{
	public:
		Radiation_gcss(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_gcss();
        void init();
        void create(Thermo<TF>&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&);

        bool check_field_exists(std::string name);
        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&);

        void exec_stats(Stats<TF>&, Thermo<TF>&, Timeloop<TF>&);
        void exec_cross(Cross<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&);
        void exec_dump(Dump<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&);
        void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&);
	private:
		//cross sections
		std::vector<std::string> crosslist;        ///< List with all crosses from ini file
		bool swcross_rflx;
        std::vector<std::string> dumplist;         ///< List with all 3d dumps from the ini file.

		void create_stats(Stats<TF>&);   ///< Initialization of the statistics.
        void create_column(Column<TF>&); ///< Initialization of the single column output.
        void create_dump(Dump<TF>&);     ///< Initialization of the single column output.
        void create_cross(Cross<TF>&);   ///< Initialization of the single column output.
        std::vector<std::string> available_masks;   // Vector with the masks that fields can provide

		using Radiation<TF>::swradiation;
		using Radiation<TF>::master;
		using Radiation<TF>::grid;
		using Radiation<TF>::fields;
		using Radiation<TF>::field3d_operators;

        TF lat;
        TF lon;
        TF xka;
        TF fr0;
        TF fr1;
        TF div;

        const TF mu_min = 0.035;


        #ifdef USECUDA
        // GPU functions and variables
        void get_radiation_field_g(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&);
        #endif

};
#endif
