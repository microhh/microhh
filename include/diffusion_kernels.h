/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef DIFFUSION_KERNELS_H
#define DIFFUSION_KERNELS_H

#include "fast_math.h"

namespace Diffusion_kernels
{
    namespace fm = Fast_math;

    template <typename TF, bool surface_model>
    void calc_strain2(
            TF* const restrict strain2,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict ugradbot,
            const TF* const restrict vgradbot,
            const TF* const restrict z,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxi, const TF dyi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;
        const int k_offset = surface_model ? 1 : 0;

        const TF zsl = z[kstart];

        // If the wall isn't resolved, calculate du/dz and dv/dz at lowest grid height using MO
        if (surface_model)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    strain2[ijk] = TF(2.)*(
                            // du/dx + du/dx
                            + fm::pow2((u[ijk+ii]-u[ijk])*dxi)

                            // dv/dy + dv/dy
                            + fm::pow2((v[ijk+jj]-v[ijk])*dyi)

                            // dw/dz + dw/dz
                            + fm::pow2((w[ijk+kk]-w[ijk])*dzi[kstart])

                            // du/dy + dv/dx
                            + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi)
                            + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi)
                            + TF(0.125)*fm::pow2((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi)
                            + TF(0.125)*fm::pow2((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi)

                            // du/dz
                            + TF(0.5) * fm::pow2(ugradbot[ij])

                            // dw/dx
                            + TF(0.125)*fm::pow2((w[ijk      ]-w[ijk-ii   ])*dxi)
                            + TF(0.125)*fm::pow2((w[ijk+ii   ]-w[ijk      ])*dxi)
                            + TF(0.125)*fm::pow2((w[ijk   +kk]-w[ijk-ii+kk])*dxi)
                            + TF(0.125)*fm::pow2((w[ijk+ii+kk]-w[ijk   +kk])*dxi)

                            // dv/dz
                            + TF(0.5) * fm::pow2(vgradbot[ij])

                            // dw/dy
                            + TF(0.125)*fm::pow2((w[ijk      ]-w[ijk-jj   ])*dyi)
                            + TF(0.125)*fm::pow2((w[ijk+jj   ]-w[ijk      ])*dyi)
                            + TF(0.125)*fm::pow2((w[ijk   +kk]-w[ijk-jj+kk])*dyi)
                            + TF(0.125)*fm::pow2((w[ijk+jj+kk]-w[ijk   +kk])*dyi) );

                    // add a small number to avoid zero divisions
                    strain2[ijk] += Constants::dsmall;
                }
        }

        for (int k=kstart+k_offset; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    strain2[ijk] = TF(2.)*(
                                   // du/dx + du/dx
                                   + fm::pow2((u[ijk+ii]-u[ijk])*dxi)

                                   // dv/dy + dv/dy
                                   + fm::pow2((v[ijk+jj]-v[ijk])*dyi)

                                   // dw/dz + dw/dz
                                   + fm::pow2((w[ijk+kk]-w[ijk])*dzi[k])

                                   // du/dy + dv/dx
                                   + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi)
                                   + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi)
                                   + TF(0.125)*fm::pow2((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi)
                                   + TF(0.125)*fm::pow2((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi)

                                   // du/dz + dw/dx
                                   + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi)
                                   + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi)
                                   + TF(0.125)*fm::pow2((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi)
                                   + TF(0.125)*fm::pow2((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi)

                                   // dv/dz + dw/dy
                                   + TF(0.125)*fm::pow2((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi)
                                   + TF(0.125)*fm::pow2((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi)
                                   + TF(0.125)*fm::pow2((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi)
                                   + TF(0.125)*fm::pow2((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi) );

                    // Add a small number to avoid zero divisions.
                    strain2[ijk] += Constants::dsmall;
                }
    }

}
#endif
