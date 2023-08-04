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

#ifndef DUST_H
#define DUST_H

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

template<typename TF>
class Dust
{
    public:
        Dust(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Dust();

        void exec(Stats<TF>&);
        void create(const double);
        unsigned long get_time_limit();

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_dust;
        TF cfl_max;
        unsigned long idt_max;

        // Gravitational settling velocities, negative downward.
        std::map<std::string, TF> ws;
};
#endif
