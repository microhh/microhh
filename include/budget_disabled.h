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

#ifndef BUDGET_DISABLED
#define BUDGET_DISABLED

#include "budget.h"

class Input;
class Master;
class Grid;
class Fields;
class Thermo;
class Diff;
class Advec;
class Stats;
struct Mask;

/**
 * Derived class for disabled budget statistics
 */
class Budget_disabled : public Budget
{
    public:
        Budget_disabled(Input*, Master*, Grid*, Fields*, Thermo*, Diff*, Advec*, Stats*);
        virtual ~Budget_disabled();

        void init();
        void create();
        void exec_stats(Mask*);
};
#endif
