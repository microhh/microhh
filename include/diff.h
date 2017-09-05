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

#ifndef DIFF
#define DIFF

// forward declaration to speed up build time
class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;

enum class Diffusion_type {Disabled, Diff_2};

template <typename TF>
class Diff
{
    public:
        Diff(Master&, Grid<TF>&, Fields<TF>&, Input&);  ///< Constructor of the diffusion class
        virtual ~Diff();                                ///< Destructor of the diffusion class

        // Pure virtual functions below which have to be implemented by the derived class
        virtual Diffusion_type get_switch() = 0;
        virtual void set_values() = 0;
        virtual void exec_viscosity() = 0;
        virtual void exec() = 0;

        virtual unsigned long get_time_limit(unsigned long, double) = 0;
        virtual double get_dn(double) = 0;

        static std::shared_ptr<Diff> factory(Master&, Grid<TF>&, Fields<TF>&, Input&, const std::string); ///< Factory function for diffusion class generation.

        //#ifdef USECUDA
        //// GPU functions and variables
        //virtual void prepare_device() = 0;
        //#endif

    protected:
        Grid<TF>& grid;
        Fields<TF>& fields;
        Master& master;
};
#endif
