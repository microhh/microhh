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

#ifndef DIFF_H
#define DIFF_H

// forward declaration to speed up build time
class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Boundary;
template<typename> class Thermo;
template<typename> class Stats;

enum class Diffusion_type {Disabled, Diff_2, Diff_4, Diff_smag2, Diff_deardorff};

template <typename TF>
class Diff
{
    public:
        Diff(Master&, Grid<TF>&, Fields<TF>&, Boundary<TF>&, Input&); ///< Constructor of the diffusion class
        virtual ~Diff(); ///< Destructor of the diffusion class

        // Pure virtual functions below which have to be implemented by the derived class
        virtual Diffusion_type get_switch() const = 0;
        virtual void create(Stats<TF>&) = 0;
        virtual void exec_viscosity(Stats<TF>&, Thermo<TF>&) = 0;
        virtual void init() = 0;
        virtual void exec(Stats<TF>&) = 0;
        virtual void exec_stats(Stats<TF>&, Thermo<TF>&) = 0;
        virtual void diff_flux(Field3d<TF>&, const Field3d<TF>&) = 0;

        virtual unsigned long get_time_limit(unsigned long, double) = 0;
        virtual double get_dn(double) = 0;

        static std::shared_ptr<Diff> factory(Master&, Grid<TF>&, Fields<TF>&, Boundary<TF>&, Input&);

        #ifdef USECUDA
        // GPU functions and variables
        virtual void prepare_device(Boundary<TF>&) = 0;
        #endif

        TF tPr;

    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Boundary<TF>& boundary;
};
#endif
