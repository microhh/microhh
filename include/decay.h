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

#ifndef DECAY_H
#define DECAY_H

#include <vector>
#include <string>
#include <map>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

/**
 * Class that creates a decay term for scalars.
 */

enum class Decay_type {disabled, enabled, exponential};

template<typename TF>
class Decay
{
    public:
        Decay(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the decay class.
        ~Decay();                                       ///< Destructor of the decay class.

        void init(Input&);           ///< Initialize the arrays that contain the profiles.
        void create(Input&, Stats<TF>&);   ///< Read the profiles of the forces from the input.
        void exec(double, Stats<TF>&);     ///< Add the tendencies belonging to the decay processes.

        void get_mask(Stats<TF>&, std::string);
        bool has_mask(std::string);


    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        // Internal switches for various forcings
        struct Decay_var
        {
            double timescale; ///< Decay timescale.
            Decay_type type; ///< Switch for the decay.
        };

        typedef std::map<std::string, Decay_var> Decay_map;
        Decay_map dmap;

        std::vector<std::string> available_masks = {"couvreux"};   // Vector with the masks that fields can provide
        TF nstd_couvreux;

        const std::string tend_name = "decay";
        const std::string tend_longname = "Decay";

};
#endif
