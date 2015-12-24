/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
class Stats;
class Grid;
class Fields;
class Thermo;
class Diff;
class Advec;
struct Mask;

/**
 * Base class for the budget statistics. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different budget statistics.
 */
class Budget
{
    public:
        Budget(Input*, Master*, Grid*, Fields*, Thermo*, Diff*, Advec*, Stats*);
        virtual ~Budget();

        static Budget* factory(Input*, Master*, Grid*, Fields*, Thermo*, Diff*, Advec*, Stats*); ///< Factory function for budget class generation.

        virtual void init() = 0;
        virtual void create() = 0;
        virtual void exec_stats(Mask*) = 0;

    protected:
        Master& master;
        Grid&   grid;
        Fields& fields;
        Thermo& thermo;
        Diff&   diff;
        Advec&  advec;
        Stats&  stats;

        std::string swbudget;
};
#endif

