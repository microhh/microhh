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

#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <vector>
#include <string>
#include <map>
#include "timeloop.h"


class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

/**
 * Class that creates a chemistry term for scalars.
 */

enum class Chemistry_type {disabled, enabled, simple};

template<typename TF>
class Chemistry
{
    public:
        Chemistry(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the chemistry class.
        ~Chemistry();                                       ///< Destructor  of the chemistry class.

        void init(Input&);                 ///< Initialize the arrays that contain the profiles.
        void create(const Timeloop<TF>&, std::string, Netcdf_handle&, Stats<TF>&);   ///< Read the profiles of the forces from the input.
        void update_time_dependent(Timeloop<TF>&); ///< Update the time dependent parameters.
        void exec(Thermo<TF>&,double,double);     ///< Add the tendencies belonging to the chemistry processes.
	void exec_stats(const int, const double, Stats<TF>&);   /// calculate statistics

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

	// internal variable
	struct Chemistry_var
	{
	     Chemistry_type type;
	};

        typedef std::map<std::string, Chemistry_var> Chemistry_map;
        Chemistry_map cmap;

	Mask<TF> m;     // borrow from Stats to gather statistics chemistry
        int statistics_counter;
	TF switch_dt;   // to switch between complicated and detailed solver
	std::vector<std::string> jname={"jo31d","jh2o2","jno2","jno3a","jno3b","jch2or","jch2om","jch3o2h"};
	std::vector<std::string> ename={"emi_isop"};
        TF jval[8];   // time-interpolated value to pass to the chemistry routine
        TF emval[1];
	std::vector<TF> time;
	std::vector<TF> jo31d;
	std::vector<TF> jh2o2;
	std::vector<TF> jno2;
	std::vector<TF> jno3a;
	std::vector<TF> jno3b;
	std::vector<TF> jch2or;
	std::vector<TF> jch2om;
	std::vector<TF> jch3o2h;
	std::vector<TF> emi_isop;
	std::vector<TF> rfa;
	TF trfa;

        const std::string tend_name = "chemistry";
        const std::string tend_longname = "Chemistry";
};
#endif
