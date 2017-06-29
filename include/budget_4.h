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

#ifndef BUDGET_4
#define BUDGET_4

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

class Budget_4 : public Budget
{
    public:
        Budget_4(Input*, Master*, Grid*, Fields*, Thermo*, Diff*, Advec*, Force*, Stats*);
        ~Budget_4();

        void init();
        void create();
        void exec_stats(Mask*);

    private:
        double* umodel;
        double* vmodel;

        void calc_ke(double*, double*, double*,
                     double*, double*,
                     double, double,
                     double*, double*);

        void calc_tke_budget_shear_turb(double*, double*, double*,
                                        double*, double*,
                                        double*, double*,
                                        double*, double*, double*, double*,
                                        double*, double*, double*, double*, double*,
                                        double*, double*);

        void calc_tke_budget(double*, double*, double*, double*,
                             double*, double*,
                             double*, double*,
                             double*, double*, double*, double*, double*,
                             double*, double*, double*, double*, double*,
                             double*, double*, double*,
                             double*, double*, double*, double*,
                             double*, double*, double);

        void calc_tke_budget_buoy(double*, double*, double*,
                                  double*, double*,
                                  double*, double*, double*);

        void calc_b2_budget(double*, double*,
                            double*,
                            double*, double*, double*, double*,
                            double*, double*,
                            double);

        void calc_bw_budget(double*, double*, double*, double*,
                            double*, double*,
                            double*, double*, double*,
                            double*, double*, double*, double*,
                            double*, double*,
                            double);

        void calc_pe(double*, double*, double*, double*,
                     double*,
                     double*,
                     double*, double*, double*,
                     double*);

        void calc_pe_budget(double*, double*, double*, double*,
                            double*, double*, double*,
                            double*, double*, double*, double*,
                            double);

        void calc_bpe_budget(double*, double*, double*, double*, double*,
                             double*, double*, double*,
                             double*,
                             double*, double*, double*,
                             double);

        double calc_zsort   (double, double*, double*, int);
        double calc_dzstardb(double, double*, double*);
};
#endif
