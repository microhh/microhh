/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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

#ifndef THERMO
#define THERMO

class Master;
class Input;
class Grid;
class Fields;
struct Mask;

/**
 * Base class for the thermo scheme. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different thermodynamics.
 */
class Thermo
{
    public:
        Thermo(Model*, Input*);
        virtual ~Thermo();
        static Thermo* factory(Master*, Input*, Model*); ///< Factory function for thermo class generation.
        std::string get_switch();

        // Below are the functions that the derived class has to implement.
        virtual void init() = 0;
        virtual void create(Input*) = 0;
        virtual void exec() = 0;

        virtual void exec_stats(Mask*) = 0;
        virtual void exec_cross() = 0;
        virtual void exec_dump() = 0;

        virtual void get_mask(Field3d*, Field3d*, Mask*) = 0;

        // Interfacing functions to get buoyancy properties from other classes.
        virtual bool check_field_exists(std::string name) = 0;
        virtual void get_thermo_field(Field3d*, Field3d*, std::string name) = 0;
        virtual void get_buoyancy_surf(Field3d*) = 0;
        virtual void get_buoyancy_fluxbot(Field3d*) = 0;
        virtual void get_prog_vars(std::vector<std::string>*) = 0;

#ifdef USECUDA
        // GPU functions and variables.
        virtual void prepare_device() = 0;
        virtual void clear_device() = 0;
#endif

    protected:
        Grid*   grid;
        Fields* fields;
        Master* master;
        Model*  model;

        std::string swthermo;
};
#endif
