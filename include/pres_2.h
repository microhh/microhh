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
#include "defines.h"

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

        std::vector<TF> bmati;
        std::vector<TF> bmatj;
        std::vector<TF> a;
        std::vector<TF> c;
        std::vector<TF> work2d;

#ifdef USECUDA
        double* bmati_g;
        double* bmatj_g;
        double* a_g;
        double* c_g;
        double* work2d_g;
#endif

        void input(double* const restrict, 
                   const double* const restrict, const double* const restrict, const double* const restrict,
                   const double* const restrict, const double* const restrict, const double* const restrict,
                   const double* const restrict, const double* const restrict, const double* const restrict,
                   const double);

        void solve(double* const restrict, double* const restrict, double*,
                   const double* const restrict, const double* const restrict,
                   double* const restrict, double* const restrict, double*, double*);

        void output(double* const restrict, double* const restrict, double* const restrict,
                    const double* const restrict, const double* const restrict);

        void tdma(double*, double*, double*, double*, 
                  double*, double*);

        double calc_divergence(const double* const restrict, const double* const restrict, const double* const restrict,
                               const double* const restrict,
                               const double* const restrict, const double* const restrict);
};
#endif
