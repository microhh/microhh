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

#ifndef MICROPHYS_DISABLED
#define MICROPHYS_DISABLED

#include "microphys.h"

class Master;
class Input;
class Data_block;

template<typename> class Grid;
template<typename> class Stats;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Thermo;
template<typename> class Field3d;
template<typename> class Microphys;

template<typename TF>
class Microphys_disabled : public Microphys<TF>
{
    public:
        Microphys_disabled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Microphys_disabled();

        void init() {};
        void create(Input&, Data_block&, Stats<TF>&, Cross<TF>&, Dump<TF>&) {};
        void exec(Thermo<TF>&, const double) {};
        void exec_stats(Stats<TF>&, std::string, Field3d<TF>&, Field3d<TF>&, const double) {};
        virtual void exec_dump(Dump<TF>&, unsigned long) {};
        virtual void exec_cross(Cross<TF>&, unsigned long) {};
        void get_mask(Field3d<TF>&, Field3d<TF>&, Stats<TF>&, std::string) {};
        bool has_mask(std::string) {return false;};

        unsigned long get_time_limit(unsigned long, double);

        #ifdef USECUDA
        //void prepare_device() {};
        //void clear_device() {};
        //void forward_device() {};
        //void backward_device() {};
        #endif

    private:
        using Microphys<TF>::swmicrophys;
};
#endif
