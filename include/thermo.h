/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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

// forward declarations to speed up build time
class cmaster;
class cgrid;
class cfields;

class cthermo
{
  public:
    cthermo(cmodel *);
    virtual ~cthermo();
    virtual int readinifile(cinput *);
    virtual int exec();
    virtual int create();

    // interfacint functions to get buoyancy properties from other classes
    virtual int getbuoyancysurf(cfield3d *);
    virtual int getbuoyancyfluxbot(cfield3d *);
    //virtual int getbuoyancy(cfield3d *, cfield3d *);
    virtual int checkthermofield(std::string name);
    virtual int getthermofield(cfield3d *, cfield3d *, std::string name);

    std::string getname();

  protected:
    cgrid   *grid;
    cfields *fields;
    cmaster *master;

    std::string swthermo;
};
#endif
