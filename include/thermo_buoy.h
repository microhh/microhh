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

#ifndef THERMO_BUOY
#define THERMO_BUOY

#include "thermo.h"

class Master;
class Grid;
class Fields;

/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */
class ThermoBuoy : public Thermo
{
  public:
    ThermoBuoy(Model *, Input *); ///< Constructor of the dry thermodynamics class.
    virtual ~ThermoBuoy();                ///< Destructor of the dry thermodynamics class.

    virtual void exec();                ///< Add the tendencies belonging to the buoyancy.

    virtual bool checkThermoField(std::string name);
    virtual void getBuoyancy(Field3d *, Field3d *);       ///< Compute the buoyancy for usage in another routine.
    virtual void getBuoyancySurf(Field3d *);              ///< Compute the near-surface and bottom buoyancy for usage in another routine.
    virtual void getBuoyancyFluxbot(Field3d *);           ///< Compute the bottom buoyancy flux for usage in another routine.
    virtual void getProgVars(std::vector<std::string> *); ///< Retrieve a list of prognostic variables.
    virtual void getThermoField(Field3d *, Field3d *, std::string name);

  private:
    void calcBuoyancy(double *, double *);         ///< Calculation of the buoyancy.
    void calcBuoyancyBot(double *, double *,
                         double *, double *);      ///< Calculation of the near-surface and surface buoyancy.
    void calcBuoyancyFluxbot(double *, double *);  ///< Calculation of the buoyancy flux at the bottom.
    void calcBuoyancyTend_2nd(double *, double *); ///< Calculation of the buoyancy tendency with 2nd order accuracy.
    void calcBuoyancyTend_4th(double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.
};
#endif
