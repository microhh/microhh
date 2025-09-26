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

        void exec();                     // Update particle velocity/tendency.
        void integrate(Timeloop<TF>&);   // Integrate particle locations.

        void create(Timeloop<TF>&);

        // Load/save restart files.
        void load(const std::string&, const int);
        void save(const std::string&, const int);

        // Dump particles.
        bool do_dump(const unsigned long);
        void dump(const int);

        unsigned long get_time_limit(const unsigned long);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_particle;       // Lagrangian particle on/off.
        int n_particles;        // Global number of particles.

        // Raw dump of all particles.
        bool sw_dump;
        unsigned long isampletime_dump;

        // Particle property arrays are oversized by a factor `buffer_margin`.
        // This reduces the amount of time that the arrays have to be resized
        // when particles move between cores.
        const TF buffer_margin = 1.2;

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
