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

#ifndef PARTICLE_LAGRANGIAN_H
#define PARTICLE_LAGRANGIAN_H

#include <hdf5.h>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Timeloop;

enum class Grid_location { Location_u, Location_v, Location_w };

template<typename TF>
class Particle_lagrangian
{
    public:
        Particle_lagrangian(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Particle_lagrangian();

        void create(Timeloop<TF>&, Stats<TF>&);    // Setup particles, I/O, statistics, ...
        void exec(Timeloop<TF>&);                  // Update particle velocity/tendency.
        void integrate(Timeloop<TF>&);             // Integrate particle locations.

        void load(const int);                      // Load particles from restart file.
        void save(const int);                      // Save particles to restart file.

        // Dump particles.
        bool do_dump(const unsigned long);         // To dump or not.
        void dump(Timeloop<TF>&);                  // Save particle dumps.

        // Statistics.
        void exec_stats(Stats<TF>&);               // Save statistics.

        unsigned long get_time_limit(const unsigned long);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_particle;                          // Lagrangian particle on/off.
        int n_particles;                           // Global number of particles.

        // Raw dump of all particles.
        bool sw_dump;
        unsigned long isampletime_dump;
        hid_t dump_file_id;                        // HDF5 file handle for particle dumps.
        int szip_compression;                      // SZIP compression pixels-per-block (0 = no compression).

        // Profile and time serie statistics.
        bool sw_stats;

        // Particle property arrays are oversized by a factor `reserve_ratio`.
        // This reduces the number of time that the arrays have to be resized
        // when particles move between cores.
        TF reserve_ratio;

        // Memory buffer increase/decrease count for statistics...
        int mem_inc = 0;
        int mem_dec = 0;

        // Delay start of particle release.
        unsigned long istarttime;

        // Particle properties.
        std::vector<int> uid;

        // Location.
        std::vector<TF> xp;
        std::vector<TF> yp;
        std::vector<TF> zp;

        // Velocity.
        std::vector<TF> up;
        std::vector<TF> vp;
        std::vector<TF> wp;

        // Tendency.
        std::vector<TF> xpt;
        std::vector<TF> ypt;
        std::vector<TF> zpt;

        // Interpolation indexes and factors.
        // Re-used for different locations on grid (u, v, w, ..).
        std::vector<int> il;
        std::vector<int> jl;
        std::vector<int> kl;

        std::vector<TF> fx;
        std::vector<TF> fy;
        std::vector<TF> fz;
    };
#endif
