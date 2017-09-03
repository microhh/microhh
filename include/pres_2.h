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

#ifndef PRES_2
#define PRES_2

#include "pres.h"

class Master;
template<typename> class Grid;
template<typename> class Fields;

template<typename TF>
class Pres_2 : public Pres<TF>
{
    public:
        Pres_2(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Pres_2();

        void init();
        void set_values();

        void exec(double);
        double check_divergence();

#ifdef USECUDA
        void prepare_device();
        void clear_device();
#endif

    private:
        using Pres<TF>::master;
        using Pres<TF>::grid;
        using Pres<TF>::fields;

        double* bmati;
        double* bmatj;
        double* a;
        double* c;
        double* work2d;

#ifdef USECUDA
        double* bmati_g;
        double* bmatj_g;
        double* a_g;
        double* c_g;
        double* work2d_g;
#endif

        void input(double*, 
                   double*, double*, double*,
                   double*, double*, double*,
                   double*, double*, double*,
                   double);

        void solve(double*, double*, double*,
                   double*, double*,
                   double*, double*, double*, double*);

        void output(double*, double*, double*,
                    double*, double*);

        void tdma(double*, double*, double*, double*, 
                  double*, double*);

        double calc_divergence(double*, double*, double*, double*, double*, double*);
};
#endif
