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

#include <algorithm>
#include <cmath>
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "thermo.h"
#include "boundary.h"
#include "stats.h"
#include "monin_obukhov.h"
#include "tools.h"

#include "diff_deardorff.h"
#include "diff_kernels.cuh"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace dk = Diff_kernels_g;

    template<typename TF, Surface_model surface_model, bool sw_mason> __global__
    void calc_evisc_g(
            TF* const restrict evisc,
            const TF* const restrict a,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict dzi,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn, const TF cm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            evisc[ijk] = 10;
        }

    //   if (surface_model == Surface_model::Disabled)
    //        throw std::runtime_error("Resolved wall not supported in Deardorff SGSm.");
    //   else
    //   {
    //        // Variables for the wall damping and length scales
    //        const TF n_mason = 2.;
    //        TF mlen;
    //        TF fac;

    //        // Calculate geometric filter width, based on Deardorff (1980)
    //        const TF mlen0 = std::pow(dx*dy*dz[kstart], TF(1./3.));

    //        for (int j=jstart; j<jend; ++j)
    //            #pragma ivdep
    //            for (int i=istart; i<iend; ++i)
    //            {
    //                const int ij = i + j*jj;
    //                const int ijk = i + j*jj + kstart*kk;

    //                if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
    //                    mlen = cn * std::sqrt(a[ijk]) / std::sqrt(bgradbot[ij]);
    //                else
    //                    mlen = mlen0;

    //                fac  = std::min(mlen0, mlen);

    //                if (sw_mason) // Apply Mason's wall correction here
    //                    fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
    //                                (std::pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);

    //                // Calculate eddy diffusivity for momentum and enforce minimum value (mvisc), as in DALES.
    //                const TF kvisc = cm * fac * std::sqrt(a[ijk]);
    //                evisc[ijk] = std::max(kvisc, mvisc<TF>);
    //            }

    //        for (int k=kstart+1; k<kend; ++k) // Counter starts at kstart (as sgstke is defined here)
    //        {
    //            // Calculate geometric filter width, based on Deardorff (1980)
    //            const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

    //            for (int j=jstart; j<jend; ++j)
    //                #pragma ivdep
    //                for (int i=istart; i<iend; ++i)
    //                {
    //                    const int ij = i + j*jj;
    //                    const int ijk = i + j*jj + k*kk;

    //                    if (N2[ijk] > 0) // Only if stably stratified, adapt length scale
    //                        mlen = cn * std::sqrt(a[ijk]) / std::sqrt(N2[ijk]);
    //                    else
    //                        mlen = mlen0;

    //                    fac  = std::min(mlen0, mlen);

    //                    if (sw_mason) // Apply Mason's wall correction here
    //                        fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
    //                                    (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

    //                    // Calculate eddy diffusivity for momentum and enforce minimum value (mvisc), as in DALES.
    //                    const TF kvisc = cm * fac * std::sqrt(a[ijk]);
    //                    evisc[ijk] = std::max(kvisc, mvisc<TF>);
    //                }
    //        }
    //    }
    }

    template<typename TF, Surface_model surface_model, bool sw_mason> __global__
    void calc_evisc_heat_g(
            TF* const restrict evisch,
            const TF* const restrict evisc,
            const TF* const restrict a,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF dx, const TF dy,
            const TF cn, const TF ch1, const TF ch2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            evisch[ijk] = 10;
        }

        //if (surface_model == Surface_model::Disabled)
        //     throw std::runtime_error("Resolved wall not supported in Deardorff SGSm.");
        //else
        //{
        //    // Variables for the wall damping and length scales
        //    const TF n_mason = 2.;
        //    TF mlen;
        //    TF fac;

        //    // Calculate geometric filter width, based on Deardorff (1980)
        //    const TF mlen0 = std::pow(dx*dy*dz[kstart], TF(1./3.));

        //    for (int j=jstart; j<jend; ++j)
        //        #pragma ivdep
        //        for (int i=istart; i<iend; ++i)
        //        {
        //            const int ij = i + j*jj;
        //            const int ijk = i + j*jj + kstart*kk;

        //            if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
        //                mlen = cn * std::sqrt(a[ijk]) / std::sqrt(bgradbot[ij]);
        //            else
        //                mlen = mlen0;

        //            fac  = std::min(mlen0, mlen);

        //            if (sw_mason) // Apply Mason's wall correction here
        //                fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
        //                            (std::pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);

        //            // Calculate eddy diffusivity for momentum and enforce minimum value (mvisc), as in DALES.
        //            const TF kvisc = (ch1 + ch2 * fac / mlen0 ) * evisc[ijk];
        //            evisch[ijk] = std::max(kvisc, mvisc<TF>);
        //        }

        //    for (int k=kstart+1; k<kend; ++k) // Counter starts at kstart (as sgstke is defined here)
        //    {
        //        // Calculate geometric filter width, based on Deardorff (1980)
        //        const TF mlen0 = std::pow(dx*dy*dz[k], TF(1./3.));

        //        for (int j=jstart; j<jend; ++j)
        //            #pragma ivdep
        //            for (int i=istart; i<iend; ++i)
        //            {
        //                const int ij = i + j*jj;
        //                const int ijk = i + j*jj + k*kk;

        //                if ( N2[ijk] > 0 ) // Only if stably stratified, adapt length scale
        //                    mlen = cn * std::sqrt(a[ijk]) / std::sqrt(N2[ijk]);
        //                else
        //                    mlen = mlen0;

        //                fac  = std::min(mlen0, mlen);

        //                if (sw_mason) // Apply Mason's wall correction here
        //                    fac = std::pow(TF(1.)/(TF(1.)/std::pow(fac, n_mason) + TF(1.)/
        //                                (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

        //                // Calculate eddy diffusivity for momentum and enforce minimum value (mvisc), as in DALES.
        //                const TF kvisc = (ch1 + ch2 * fac / mlen0 ) * evisc[ijk];
        //                evisch[ijk] = std::max(kvisc, mvisc<TF>);
        //            }
        //    }
        //}
    }


}

#ifdef USECUDA
template<typename TF>
unsigned long Diff_deardorff<TF>::get_time_limit(const unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
double Diff_deardorff<TF>::get_dn(const double dt)
{
    return -1;
}

template<typename TF>
void Diff_deardorff<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = TF(1)/(gd.dx * gd.dx);
    const TF dyidyi = TF(1)/(gd.dy * gd.dy);

    // Dummy tPr value for `diff_c`.
    const TF tPr_i_dummy = 1;

    dk::diff_uvw_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
            fields.mt.at("u")->fld_g,
            fields.mt.at("v")->fld_g,
            fields.mt.at("w")->fld_g,
            fields.sd.at("evisc")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            fields.mp.at("u")->flux_bot_g,
            fields.mp.at("u")->flux_top_g,
            fields.mp.at("v")->flux_bot_g,
            fields.mp.at("v")->flux_top_g,
            gd.dzi_g,
            gd.dzhi_g,
            fields.rhoref_g,
            fields.rhorefh_g,
            gd.dxi,
            gd.dyi,
            fields.visc,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    for (auto it : fields.st)
    {
        TF* evisc_ptr;

        if (it.first == "sgstke")  // sgstke diffuses with eddy viscosity for momentum
            evisc_ptr = fields.sd.at("evisc")->fld_g;
        else  // all other scalars, normally diffuse with eddy viscosity for heat/scalars
        {
            if (!sw_buoy) // but not if there is no buoyancy (then eviscs not defined)
                evisc_ptr = fields.sd.at("evisc")->fld_g;
            else
                evisc_ptr = fields.sd.at("eviscs")->fld_g;
        }

        dk::diff_c_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
                it.second->fld_g,
                fields.sp.at(it.first)->fld_g,
                evisc_ptr,
                fields.sp.at(it.first)->flux_bot_g,
                fields.sp.at(it.first)->flux_top_g,
                gd.dzi_g,
                gd.dzhi_g,
                fields.rhoref_g,
                fields.rhorefh_g,
                dxidxi,
                dyidyi,
                tPr_i_dummy,
                fields.sp.at(it.first)->visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
}

template<typename TF>
void Diff_deardorff<TF>::exec_viscosity(Stats<TF>& stats, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    // Contain the full icells and jcells in this grid.
    const int grid2di  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int grid2dj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 grid2dGPU (grid2di, grid2dj);
    dim3 block2dGPU(blocki, blockj);

    // Get MO gradients velocity and roughness length:
    TF* dudz_g  = boundary.get_dudz_g();
    TF* dvdz_g  = boundary.get_dvdz_g();
    TF* z0m_g   = boundary.get_z0m_g();

    auto str2_tmp = fields.get_tmp_g();

    // Calculate total strain rate
    dk::calc_strain2_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
        str2_tmp->fld_g,
        fields.mp.at("u")->fld_g,
        fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g,
        dudz_g, dvdz_g,
        gd.dzi_g, gd.dzhi_g,
        gd.dxi, gd.dyi,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Start with retrieving the stability information
    if (!sw_buoy)
    {
    //    auto evisc_neutral_wrapper = [&]<Surface_model surface_model, bool sw_mason>()
    //    {
    //        calc_evisc_neutral<TF, surface_model, sw_mason>(
    //                fields.sd.at("evisc")->fld.data(),
    //                fields.sp.at("sgstke")->fld.data(),
    //                fields.mp.at("u")->fld.data(),
    //                fields.mp.at("v")->fld.data(),
    //                fields.mp.at("w")->fld.data(),
    //                gd.z.data(), gd.dz.data(),
    //                gd.dzhi.data(), z0m.data(),
    //                gd.dx, gd.dy, gd.zsize,
    //                this->cm, this->cn, fields.visc,
    //                gd.istart, gd.iend,
    //                gd.jstart, gd.jend,
    //                gd.kstart, gd.kend,
    //                gd.icells, gd.jcells,
    //                gd.ijcells,
    //                boundary_cyclic);
    //    };

    //    // Calculate eddy viscosity using MO at lowest model level
    //    if (sw_mason)
    //        evisc_neutral_wrapper.template operator()<Surface_model::Enabled, true>();
    //    else
    //        evisc_neutral_wrapper.template operator()<Surface_model::Enabled, false>();

    //    sgstke_diss_tend_neutral<TF>(
    //            fields.st.at("sgstke")->fld.data(),
    //            fields.sp.at("sgstke")->fld.data(),
    //            gd.z.data(),
    //            gd.dz.data(),
    //            z0m.data(),
    //            gd.dx,
    //            gd.dy,
    //            this->ce1,
    //            this->ce2,
    //            gd.istart, gd.iend,
    //            gd.jstart, gd.jend,
    //            gd.kstart, gd.kend,
    //            gd.icells, gd.ijcells,
    //            sw_mason);

    //            stats.calc_tend(*fields.st.at("sgstke"), tend_name_diss);
    }
    else
    {
        // Assume buoyancy calculation is needed
        auto buoy_tmp = fields.get_tmp_g();
        thermo.get_thermo_field_g(*buoy_tmp, "N2", false);

        // Get MO gradient buoyancy:
        TF* dbdz_g  = boundary.get_dbdz_g();

        // Note BvS: templated lambda functions are not (yet?) allowed by NVCC :-(
        if (sw_mason)
            calc_evisc_g<TF, Surface_model::Enabled, true><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("v")->fld_g,
                    fields.mp.at("w")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g, gd.dzi_g, z0m_g,
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells,
                    gd.ijcells);
        else
            calc_evisc_g<TF, Surface_model::Enabled, false><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("v")->fld_g,
                    fields.mp.at("w")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g, gd.dzi_g, z0m_g,
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells,
                    gd.ijcells);

        if (sw_mason)
            calc_evisc_heat_g<TF, Surface_model::Enabled, true><<<gridGPU, blockGPU>>>(
                    fields.sd.at("eviscs")->fld_g,
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g, z0m_g,
                    gd.dx, gd.dy,
                    this->cn, this->ch1, this->ch2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells);
        else
            calc_evisc_heat_g<TF, Surface_model::Enabled, false><<<gridGPU, blockGPU>>>(
                    fields.sd.at("eviscs")->fld_g,
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g, z0m_g,
                    gd.dx, gd.dy,
                    this->cn, this->ch1, this->ch2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells);

        boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);
        boundary_cyclic.exec_g(fields.sd.at("eviscs")->fld_g);

        //// BvS: I left the tendency calculations of sgstke here; feels a bit strange
        //// to calculate them in `exec_viscosity`, but otherwise strain^2 has to be
        //// recalculated in diff->exec()...
        //sgstke_buoy_tend<TF>(
        //        fields.st.at("sgstke")->fld.data(),
        //        fields.sp.at("sgstke")->fld.data(),
        //        fields.sd.at("eviscs")->fld.data(),
        //        buoy_tmp->fld.data(),
        //        dbdz.data(),
        //        gd.istart, gd.iend,
        //        gd.jstart, gd.jend,
        //        gd.kstart, gd.kend,
        //        gd.icells, gd.ijcells);

        //stats.calc_tend(*fields.st.at("sgstke"), tend_name_buoy);

        //sgstke_diss_tend<TF>(
        //        fields.st.at("sgstke")->fld.data(),
        //        fields.sp.at("sgstke")->fld.data(),
        //        buoy_tmp->fld.data(),
        //        dbdz.data(),
        //        gd.z.data(),
        //        gd.dz.data(),
        //        z0m.data(),
        //        gd.dx, gd.dy,
        //        this->cn,
        //        this->ce1,
        //        this->ce2,
        //        gd.istart, gd.iend,
        //        gd.jstart, gd.jend,
        //        gd.kstart, gd.kend,
        //        gd.icells, gd.ijcells, sw_mason); ///< SvdL, 14-04-2023: not the nicest, see above);

        //stats.calc_tend(*fields.st.at("sgstke"), tend_name_diss);

        fields.release_tmp_g(buoy_tmp);
    }

    //sgstke_shear_tend<TF>(
    //        fields.st.at("sgstke")->fld.data(),
    //        fields.sp.at("sgstke")->fld.data(),
    //        fields.sd.at("evisc")->fld.data(),
    //        str2_tmp->fld.data(),
    //        gd.istart, gd.iend,
    //        gd.jstart, gd.jend,
    //        gd.kstart, gd.kend,
    //        gd.icells, gd.ijcells);

    //stats.calc_tend(*fields.st.at("sgstke"), tend_name_shear);

    // Release temporary fields
    fields.release_tmp_g(str2_tmp);
}

template<typename TF>
void Diff_deardorff<TF>::prepare_device(Boundary<TF>& boundary)
{
}

template<typename TF>
void Diff_deardorff<TF>::clear_device()
{
}
#endif

template class Diff_deardorff<double>;
template class Diff_deardorff<float>;
