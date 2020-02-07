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

#ifndef SOIL
#define SOIL

class Master;
class Input;

enum class Soil_type {Disabled, Fixed, Interactive};

/**
 * Base class for the soil scheme. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different microphysiscs schemes.
 */
template<typename TF>
class Soil
{
    public:
        Soil(Master&, Input&);
        virtual ~Soil();

        static std::shared_ptr<Soil> factory(Master&, Input&);
        Soil_type get_switch();

        // Below are the functions that the derived class has to implement.
        virtual void init() = 0;
        virtual void create(Input&, Netcdf_handle&) = 0;
        virtual void exec() = 0;

    protected:
        Master& master;
        //Fields<TF>& fields;

        Soil_type sw_soil;
};
#endif
