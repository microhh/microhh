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

#include <cstdio>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "immersed_boundary.h"

Immersed_boundary::Immersed_boundary(Master& masterin, Grid& gridin) : 
    master(masterin),
    grid(gridin)
{
}

Immersed_boundary::~Immersed_boundary()
{
}

namespace
{
    void ib(double* const restrict u, double* const restrict w,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const int ibc_istart = istart + 7 * (iend-istart) / 64;
        const int ibc_iend   = istart + 9 * (iend-istart) / 64;
        const int ibc_kstart = kstart + 7 * (kend-kstart) / 16;
        const int ibc_kend   = kstart + 9 * (kend-kstart) / 16;

        // Set the u ghost cells
        for (int k=ibc_kstart; k<ibc_kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=ibc_istart; i<ibc_iend+1; ++i)
                {
                    const int ijk  = i + j*jj + k*kk;
                    u[ijk] = 0.;
                }

        // Set the w ghost cells
        for (int k=ibc_kstart; k<ibc_kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=ibc_istart; i<ibc_iend; ++i)
                {
                    const int ijk  = i + j*jj + k*kk;
                    w[ijk] = 0.;
                }
    }
}

void Immersed_boundary::set_ghost_cells(Fields& fields)
{
    ib(fields.u->data, fields.w->data,
       grid.istart, grid.iend,
       grid.jstart, grid.jend,
       grid.kstart, grid.kend,
       grid.icells, grid.ijcells);
}
