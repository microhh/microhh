/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef SOURCE_KERNELS_H
#define SOURCE_KERNELS_H

template<typename TF>
std::vector<int> calc_shape(
        const TF* restrict x, const TF x0, const TF sigma_x, const TF line_x, int istart, int iend)
{
    std::vector<int> range(2);

    int i = istart;
    range[0] = iend;

    for (; i<iend; ++i)
    {
        if ( x[i]-x0 + TF(4)*sigma_x > TF(0) )
        {
            range[0] = i;
            break;
        }
    }

    i = istart;
    for (; i<iend; ++i)
    {
        range[1] = iend;

        if ( x[i]-x0-line_x - TF(4)*sigma_x > TF(0) )
        {
            range[1] = i;
            break;
        }
    }

    return range;
}

#endif
