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

#ifndef DIFF
#define DIFF

// forward declaration to speed up build time
class Model;
class Grid;
class Fields;
class Master;

class Diff
{
    public:
        Diff(Model*, Input*);
        virtual ~Diff();
        static Diff* factory(Master*, Input*, Model*, const std::string); ///< Factory function for diff class generation.

        std::string get_name();

        // Pure virtual functions below.
        virtual void set_values() = 0;
        virtual void exec_viscosity() = 0;
        virtual void exec() = 0;

        virtual unsigned long get_time_limit(unsigned long, double) = 0;
        virtual double get_dn(double) = 0;

        #ifdef USECUDA
        // GPU functions and variables
        virtual void prepare_device() = 0;
        #endif

    protected:
        Model*  model;
        Grid*   grid;
        Fields* fields;
        Master* master;

        std::string swdiff;

        double dnmax;

};
#endif
