/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#ifndef BUDGET
#define BUDGET

#include <string>

class Input;
class Master;
template<typename> class Stats;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Thermo;
template<typename> class Diff;
template<typename> class Advec;
template<typename> class Force;
//struct Mask;

/**
 * Base class for the budget statistics. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different budget statistics.
 */
template<typename TF>
class Budget
{
    public:
        Budget(Master&, Grid<TF>&, Fields<TF>&, Diff<TF>&, Advec<TF>&, Force<TF>&, Stats<TF>&, Input&);
        virtual ~Budget();

        static std::shared_ptr<Budget> factory(Master&, Grid<TF>&, Fields<TF>&, Diff<TF>&, Advec<TF>&, Force<TF>&, Stats<TF>&, Input&); ///< Factory function for budget class generation.

        //virtual void init() = 0;
        //virtual void create() = 0;
        //virtual void exec_stats(Mask*) = 0;

    protected:
        Master& master;
        Grid<TF>&   grid;
        Fields<TF>& fields;
        //Thermo<TF>& thermo;
        Diff<TF>&   diff;
        Advec<TF>&  advec;
        Force<TF>&  force;
        Stats<TF>&  stats;

        enum class Budget_type {Disabled, Second_order, Fourth_order};
        Budget_type swbudget;
};
#endif

