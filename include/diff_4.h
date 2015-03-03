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

#ifndef DIFF_4
#define DIFF_4

#include "diff.h"
#include "defines.h"

class Diff4 : public Diff
{
  public:
    Diff4(Model *, Input *);
    ~Diff4();

    void setValues();
    void exec();

    unsigned long getTimeLimit(unsigned long, double);
    double get_dn(double);

  private:
    double dnmul;

    template<bool> void diffc(double * restrict, double * restrict, double * restrict, double * restrict, double);
    template<bool> void diffw(double * restrict, double * restrict, double * restrict, double * restrict, double);
};
#endif
