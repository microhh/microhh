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

#include <cstdio>
#include <cmath>
#include <stdlib.h>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "defines.h"
#include "data_block.h"


namespace
{
    template<typename TF>
    void calc_buffer(TF* const restrict at, const TF* const restrict a,
                        const TF* const restrict abuf, const TF* const restrict z, const TF zstart, const TF zsize, const TF beta, const TF sigma, const int istart, const int iend, const int icells, const int jstart, const int jend, const int ijcells, const int bufferkstart, const int kend)
    {
        const TF zsizebuf = zsize - zstart;

        TF sigmaz;

        for (int k=bufferkstart; k<kend; ++k)
        {
            sigmaz = sigma*std::pow((z[k]-zstart)/zsizebuf, beta);
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    at[ijk] -= sigmaz*(a[ijk]-abuf[k]);
                }
        }
    }

}
template<typename TF>
Buffer<TF>::Buffer(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{

    swbuffer = inputin.get_item<bool>("buffer", "swbuffer", "", false);

    if (swbuffer)
    {
        swupdate = inputin.get_item<bool>("buffer", "swupdate", "", false);
        zstart   = inputin.get_item<TF>("buffer", "zstart", "");
        sigma    = inputin.get_item<TF>("buffer", "sigma", "",2.);
        beta     = inputin.get_item<TF>("buffer", "beta", "",2.);
    }


    if (swbuffer && swupdate)
        fields.set_calc_mean_profs(true);
}

template<typename TF>
Buffer<TF>::~Buffer()
{
    #ifdef USECUDA
    clear_device();
    #endif
}

template<typename TF>
void Buffer<TF>::init()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (swbuffer)
    {
        if (swupdate)
        {
            bufferprofs["w"].resize(gd.kcells);
        }
        else
        {
            // Allocate the buffer arrays.
            for (auto& it : fields.ap )
                bufferprofs[it.first].resize(gd.kcells);
        }
    }
}

template<typename TF>
void Buffer<TF>::create(Input& inputin, Data_block& profs)
{

    if (swbuffer)
    {
        const Grid_data<TF>& gd = grid.get_grid_data();

        // Find the starting points.
        bufferkstart  = gd.kstart;
        bufferkstarth = gd.kstart;

        for (int k=gd.kstart; k<gd.kend; ++k)
        {
            // Check if the cell center is in the buffer zone.
            if (gd.z[k] < zstart)
                ++bufferkstart;
            // Check if the cell face is in the buffer zone.
            if (gd.zh[k] < zstart)
                ++bufferkstarth;
        }

        // Check whether the lowest of the two levels is contained in the buffer layer.
        if (bufferkstarth == gd.kend)
        {
            master.print_error("buffer is too close to the model top\n");
        }

        // Allocate the buffer for w on 0.
        for (int k=0; k<gd.kcells; ++k)
             bufferprofs["w"][k] = 0.;

        if (!swupdate)
        {
            // Set the buffers according to the initial profiles of the variables.
            profs.get_vector(bufferprofs["u"], "u", gd.kmax, 0, gd.kstart);
            profs.get_vector(bufferprofs["v"], "v", gd.kmax, 0, gd.kstart);
            // In case of u and v, subtract the grid velocity.
            for (int k=gd.kstart; k<gd.kend; ++k)
            {
                bufferprofs["u"][k] -= grid.utrans;
                bufferprofs["v"][k] -= grid.vtrans;
            }


            for (auto& it : fields.sp)
                profs.get_vector(bufferprofs[it.first], it.first, gd.kmax, 0, gd.kstart);
        }
    }
}

#ifndef USECUDA
template<typename TF>
void Buffer<TF>::exec()
{
    if (swbuffer)
    {
        const Grid_data<TF>& gd = grid.get_grid_data();

        if (swupdate)
        {
            // Calculate the buffer tendencies.
            calc_buffer(fields.mt.at("u")->fld.data(),fields.mp.at("u")->fld.data(),fields.mp.at("u")->fld_mean.data(),gd.z.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstart, gd.kend);
            calc_buffer(fields.mt.at("v")->fld.data(),fields.mp.at("v")->fld.data(),fields.mp.at("v")->fld_mean.data(),gd.z.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstart, gd.kend);
            calc_buffer(fields.mt.at("w")->fld.data(),fields.mp.at("w")->fld.data(),fields.mp.at("w")->fld_mean.data(),gd.zh.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstarth, gd.kend);

            for (auto& it : fields.sp)
                calc_buffer(fields.st.at(it.first)->fld.data(),fields.sp.at(it.first)->fld.data(),fields.sp.at(it.first)->fld_mean.data(),gd.z.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstart, gd.kend);
        }
        else
        {
            // Calculate the buffer tendencies.
            calc_buffer(fields.mt.at("u")->fld.data(),fields.mp.at("u")->fld.data(),bufferprofs.at("u").data(),gd.z.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstart, gd.kend);
            calc_buffer(fields.mt.at("v")->fld.data(),fields.mp.at("v")->fld.data(),bufferprofs.at("v").data(),gd.z.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstart, gd.kend);
            calc_buffer(fields.mt.at("w")->fld.data(),fields.mp.at("w")->fld.data(),bufferprofs.at("w").data(),gd.zh.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstarth, gd.kend);

            for (auto& it : fields.sp)
                calc_buffer(fields.st.at(it.first)->fld.data(),bufferprofs.at(it.first).data(),fields.sp.at(it.first)->fld_mean.data(),gd.z.data(), zstart, gd.zsize, beta, sigma, gd.istart, gd.iend, gd.icells, gd.jstart, gd.jend, gd.ijcells, bufferkstart, gd.kend);
        }
    }
}
#endif

template class Buffer<double>;
template class Buffer<float>;
