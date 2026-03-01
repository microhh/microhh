
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

#include <iostream>

#include "chemistry.h"
#include "chemistry_plume_kernels.cuh"
#include "deposition.h"

#include "tools.h"
#include "fields.h"
#include "thermo.h"
#include "constants.h"

namespace cpkg = Chemistry_plume_kernels_g;


#ifdef USECUDA
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo, const double sdt, const double dt)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    // Calculate the mean temperature profile.
    auto temperature = fields.get_tmp_g();
    thermo.get_thermo_field_g(*temperature, "T", true);
    field3d_operators.calc_mean_profile_g(temperature->fld_mean_g, temperature->fld_g);

    cpkg::pss<TF><<<gridGPU, blockGPU>>>(
        fields.st.at("hno3")->fld_g,
        fields.st.at("h2o2")->fld_g,
        fields.st.at("co")  ->fld_g,
        fields.st.at("hcho")->fld_g,
        fields.st.at("rooh")->fld_g,
        fields.st.at("c3h6")->fld_g,
        fields.st.at("o3")  ->fld_g,
        fields.st.at("no")  ->fld_g,
        fields.st.at("no2") ->fld_g,
        fields.sp.at("hno3")->fld_g,
        fields.sp.at("h2o2")->fld_g,
        fields.sp.at("co")  ->fld_g,
        fields.sp.at("hcho")->fld_g,
        fields.sp.at("rooh")->fld_g,
        fields.sp.at("c3h6")->fld_g,
        fields.sp.at("o3")  ->fld_g,
        fields.sp.at("no")  ->fld_g,
        fields.sp.at("no2") ->fld_g,
        jval_g,
        vdo3_g,
        vdno_g,
        vdno2_g,
        vdhno3_g,
        vdh2o2_g,
        vdrooh_g,
        vdhcho_g,
        temperature->fld_mean_g,
        fields.sp.at("qt")->fld_mean_g,
        gd.dzi_g,
        fields.rhoref_g,
        sdt,
        gd.istart,
        gd.iend,
        gd.jstart,
        gd.jend,
        gd.kstart,
        gd.kend,
        gd.icells,
        gd.ijcells);

    fields.release_tmp_g(temperature);
}

template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary)
{
    if (!sw_chemistry)
        return;

    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);

    // There is nothing the GPU can speedup here, simply calculate at CPU and memcpy to device.
    jval[Jval::o31d]   = ifac.fac0 * jo31d[ifac.index0]   + ifac.fac1 * jo31d[ifac.index1];
    jval[Jval::h2o2]   = ifac.fac0 * jh2o2[ifac.index0]   + ifac.fac1 * jh2o2[ifac.index1];
    jval[Jval::no2]    = ifac.fac0 * jno2[ifac.index0]    + ifac.fac1 * jno2[ifac.index1];
    jval[Jval::no3]    = ifac.fac0 * jno3[ifac.index0]    + ifac.fac1 * jno3[ifac.index1];
    jval[Jval::n2o5]   = ifac.fac0 * jn2o5[ifac.index0]   + ifac.fac1 * jn2o5[ifac.index1];
    jval[Jval::ch2or]  = ifac.fac0 * jch2or[ifac.index0]  + ifac.fac1 * jch2or[ifac.index1];
    jval[Jval::ch2om]  = ifac.fac0 * jch2om[ifac.index0]  + ifac.fac1 * jch2om[ifac.index1];
    jval[Jval::ch3o2h] = ifac.fac0 * jch3o2h[ifac.index0] + ifac.fac1 * jch3o2h[ifac.index1];

    cuda_safe_call(cudaMemcpy(jval_g, jval.data(), n_jval*sizeof(TF), cudaMemcpyHostToDevice));

    deposition->update_time_dependent(
            timeloop,
            boundary,
            vdo3_g,
            vdno_g,
            vdno2_g,
            vdhno3_g,
            vdh2o2_g,
            vdrooh_g,
            vdhcho_g);
}

template<typename TF>
void Chemistry<TF>::prepare_device()
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    const int time_size = time.size();

    // Allocate GPU arrays.
    // Only `jval_g` is enough; interpolation is performed at CPU, and memcpy'd to device.
    jval_g.allocate(n_jval);

    vdo3_g.allocate(gd.ijcells);
    vdno_g.allocate(gd.ijcells);
    vdno2_g.allocate(gd.ijcells);
    vdhno3_g.allocate(gd.ijcells);
    vdh2o2_g.allocate(gd.ijcells);
    vdrooh_g.allocate(gd.ijcells);
    vdhcho_g.allocate(gd.ijcells);

    // Copy data from host to device.
    cuda_safe_call(cudaMemcpy(jval_g,    jval.data(),    n_jval*sizeof(TF),     cudaMemcpyHostToDevice));

    const int memsize_ij = gd.ijcells * sizeof(TF);
    cuda_safe_call(cudaMemcpy(vdo3_g,    vdo3.data(),    memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdno_g,    vdno.data(),    memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdno2_g,   vdno2.data(),   memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdhno3_g,  vdhno3.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdh2o2_g,  vdh2o2.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdrooh_g,  vdrooh.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdhcho_g,  vdhcho.data(),  memsize_ij, cudaMemcpyHostToDevice));

    // Prepare deposition.
    deposition->prepare_device();
}

template<typename TF>
void Chemistry<TF>::backward_device()
{
    if (!sw_chemistry)
        return;

    deposition->backward_device();
}

template<typename TF>
void Chemistry<TF>::clear_device()
{
    if (!sw_chemistry)
        return;
}
#endif

#ifdef FLOAT_SINGLE
template class Chemistry<float>;
#else
template class Chemistry<double>;
#endif