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

#ifndef PARTICLE_BIN_H
#define PARTICLE_BIN_H

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Timeloop;

template<typename TF>
class Particle_bin
{
    public:
        Particle_bin(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Particle_bin();

        void exec(Stats<TF>&);
        void create(Timeloop<TF>&);
        unsigned long get_time_limit();

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_particle;
        TF cfl_max;
        unsigned long idt_max;

        // Gravitational settling velocities, negative downward.
        std::map<std::string, TF> w_particle;
};
#endif
