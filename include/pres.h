/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

#ifndef PRES
#define PRES

class Model;
class Grid;
class Fields;
class Master;

class Pres
{
  public:
    Pres(Model *, Input *);
    virtual ~Pres();
    static Pres* factory(Master *, Input *, Model *, const std::string); ///< Factory function for pres class generation.

    virtual void init();
    virtual void setValues();

    virtual void exec(double);
    virtual double checkDivergence();

    virtual int prepareDevice();

  protected:
    Master *master;
    Model  *model;
    Grid   *grid;
    Fields *fields;
};
#endif
