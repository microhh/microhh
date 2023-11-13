/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef BUDGET_4_H
#define BUDGET_4_H

#include "budget.h"

template<typename> class Field3d_operators;

template<typename TF>
class Budget_4 : public Budget<TF>
{
    public:
        Budget_4(Master&, Grid<TF>&, Fields<TF>&, Thermo<TF>&, Diff<TF>&, Advec<TF>&, Force<TF>&, Input&);
        ~Budget_4();

        void init();
        void create(Stats<TF>&);
        void exec_stats(Stats<TF>&);

    private:
        using Budget<TF>::master;
        using Budget<TF>::grid;
        using Budget<TF>::fields;
        using Budget<TF>::thermo;
        using Budget<TF>::diff;
        using Budget<TF>::advec;
        using Budget<TF>::force;

        Field3d_operators<TF> field3d_operators;

        std::vector<TF> umodel;
        std::vector<TF> vmodel;
        std::vector<TF> wmodel;
};
#endif
