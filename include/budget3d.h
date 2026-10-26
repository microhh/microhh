/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#ifndef BUDGET3D_H
#define BUDGET3D_H

#include <string>
#include <vector>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Thermo;
template<typename> class Dump;

template<typename TF>
class Budget3d
{
    public:
        Budget3d(Master&, Grid<TF>&, Fields<TF>&, Thermo<TF>&, Dump<TF>&, Input&);
        ~Budget3d();

        void exec(const double sub_dt, const double dt, const bool in_substep);   ///< Accumulate one RK substep's tendency.
        void prepare_dump();              ///< Normalize accumulators into a mean tendency before writing.
        void exec_dump(Dump<TF>&, unsigned long);
        void reset();                     ///< Zero the accumulators after writing.

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Thermo<TF>& thermo;

        bool swbudget3d;
        std::vector<std::string> tendencylist;
        double elapsed_time;

        std::vector<std::string> fluxlist;
        std::string evisc_name;
        TF evisc_fac;
};
#endif
