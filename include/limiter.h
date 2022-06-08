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

#ifndef LIMITER_H
#define LIMITER_H

#include <vector>
#include <string>
#include <map>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Diff; // tentativechange, SvdL, 07.06.22

template<typename TF>
class Limiter
{
    public:
        Limiter(Master&, Grid<TF>&, Fields<TF>&, Diff<TF>&, Input&); // Constructor of the decay class. // tentativechange, SvdL, 07.06.22
        ~Limiter();                                       // Destructor of the decay class.

        void create(Stats<TF>&); // Read the profiles of the forces from the input.
        void exec(double, Stats<TF>&); // Add the tendencies belonging to the decay processes.

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        std::vector<std::string> limit_list;

        const std::string tend_name = "limit";
        const std::string tend_longname = "Limiter";
};
#endif
