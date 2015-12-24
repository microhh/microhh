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

#ifndef BUDGET_2
#define BUDGET_2

class Input;
class Master;
class Stats;
class Grid;
class Fields;
class Thermo;
class Diff;
class Advec;
struct Mask;

#include "budget.h"

class Budget_2 : public Budget
{
    public:
        Budget_2(Input*, Master*, Grid*, Fields*, Thermo*, Diff*, Advec*, Stats*);
        ~Budget_2();

        void init();
        void create();
        void exec_stats(Mask*);

    private:
        double* umodel;
        double* vmodel;

        void calc_kinetic_energy(double*, double*, const double*, const double*, const double*, const double*, const double*, const double, const double);

        void calc_advection_terms(double*, double*, double*, double*, double*, double*, double*, double*, double*,
                                  const double*, const double*, const double*, const double*, const double*,
                                  double*, double*, const double*, const double*); 

        void calc_pressure_terms(double*, double*, double*, double*, double*, double*, double*,
                                 const double*, const double*, const double*, const double*,
                                 const double*, const double*, const double*, const double*,
                                 const double, const double);

        void calc_diffusion_terms_DNS(double*, double*, double*, double*, double*, double*,
                                      double*, double*, double*, double*, double*, double*, double*,
                                      const double*, const double*, const double*, const double*,
                                      const double*, const double*, const double*,
                                      const double, const double, const double);

        void calc_diffusion_terms_LES(double*, double*, double*, double*, double*, double*,
                                      double*, double*, double*, double*, double*, double*, double*,
                                      const double*, const double*, const double*, const double*, const double*,
                                      const double*, const double*, const double*, const double*, const double*,
                                      const double, const double);

        void calc_buoyancy_terms(double*, double*, double*,
                                 const double*, const double*, const double*, const double*, const double*);
};
#endif
