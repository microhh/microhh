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

class Thermo
{
  public:
    Thermo(Model *, Input *);
    virtual ~Thermo();
    static Thermo* factory(Master *, Input *, Model *); ///< Factory function for thermo class generation.

    virtual void init();
    virtual void create(Input *);
    virtual void exec();
    virtual void exec_stats(Mask *);

    virtual void exec_cross();
    virtual void exec_dump();

    virtual void get_mask(Field3d *, Field3d *, Mask *);

    // interfacing functions to get buoyancy properties from other classes
    virtual bool checkThermoField(std::string name);
    virtual void getThermoField(Field3d *, Field3d *, std::string name);
    virtual void getBuoyancySurf(Field3d *);
    virtual void getBuoyancyFluxbot(Field3d *);
    virtual void getProgVars(std::vector<std::string> *);

    std::string getSwitch();

    // GPU functions and variables
    virtual void prepare_device();
    virtual void clear_device();

  protected:
    Grid   *grid;
    Fields *fields;
    Master *master;
    Model  *model;

    std::string swthermo;
};
#endif
