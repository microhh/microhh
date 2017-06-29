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

#ifndef DIFF_2
#define DIFF_2

#include "diff.h"

class Diff_2 : public Diff
{
    public:
        Diff_2(Model*, Input*);
        ~Diff_2();

        void set_values();
        void exec();

        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        // Empty functions, these are allowed to pass.
        void exec_viscosity() {}

        #ifdef USECUDA
        void prepare_device() {};
        #endif

    private:
        double dnmul;

        void diff_c(double*, double*, double*, double*, double);
        void diff_w(double*, double*, double*, double*, double);
};
#endif
