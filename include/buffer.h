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

#ifndef BUFFER_H
#define BUFFER_H

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

/**
 * Class for the buffer layer in the top of the domain.
 * This class performs the gravity wave damping in the top of the domain to
 * prevent reflection at the top boundary.
 */
 template<typename TF>
class Buffer
{
    public:
        Buffer(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the buffer class.
        ~Buffer(); ///< Destructor of the buffer class.

        void init(); ///< Initialize the arrays that contain the profiles.
        void create(Input&, Netcdf_handle&, Stats<TF>&); ///< Read the profiles of the forces from the input.
        void exec(Stats<TF>&); ///< Add the tendencies created by the damping.

        // GPU functions and variables
        void prepare_device(); ///< Allocate and copy buffer profiles at/to GPU
        void clear_device(); ///< Allocate and copy buffer profiles at/to GPU

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        TF zstart; ///< Height above which the buffer layer is starting.
        TF sigma;  ///< Damping frequency.
        TF beta;   ///< Exponent for damping increase with height.

        int bufferkstart;  ///< Grid point at cell center at which damping starts.
        int bufferkstarth; ///< Grid point at cell face at which damping starts.

        std::map<std::string, std::vector<TF>> bufferprofs; ///< Map containing the buffer profiles.

        bool swbuffer; ///< Switch for buffer.
        bool swupdate; ///< Switch for enabling runtime updating of buffer profile.

        // GPU functions and variables
        std::map<std::string, TF*> bufferprofs_g; ///< Map containing the buffer profiles at GPU.

        const std::string tend_name = "damp";
        const std::string tend_longname = "Damping";
};
#endif
