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

#include "grid.h"
#include "fields.h"
#include "pres_4.h"

#ifdef USECUDA
void cpres_4::exec(double dt)
{
  // 1. Create the input for the pressure solver.
  pres_in(fields->sd["p"]->data,
          fields->u ->data, fields->v ->data, fields->w ->data,
          fields->ut->data, fields->vt->data, fields->wt->data, 
          grid->dzi4, dt);

  // 2. Solve the Poisson equation using FFTs and a heptadiagonal solver
  // Take slices out of a temporary field to save memory. The temp arrays
  // are always big enough, this cannot fail.
  double *tmp2 = fields->sd["tmp2"]->data;
  double *tmp3 = fields->sd["tmp3"]->data;
  const int ns = grid->iblock*(grid->kmax+4);
  pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, grid->dz,
             m1, m2, m3, m4,
             m5, m6, m7,
             &tmp2[0*ns], &tmp2[1*ns], &tmp2[2*ns], &tmp2[3*ns], 
             &tmp3[0*ns], &tmp3[1*ns], &tmp3[2*ns], &tmp3[3*ns], 
             bmati, bmatj);

  // 3. Get the pressure tendencies from the pressure field.
  pres_out(fields->ut->data, fields->vt->data, fields->wt->data, 
           fields->sd["p"]->data, grid->dzhi4);
}
#endif
