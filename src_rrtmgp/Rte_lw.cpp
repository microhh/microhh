/*
 * This file is part of a C++ interface to the Radiative Transfer for Energetics (RTE)
 * and Rapid Radiative Transfer Model for GCM applications Parallel (RRTMGP).
 *
 * The original code is found at https://github.com/RobertPincus/rte-rrtmgp.
 *
 * Contacts: Robert Pincus and Eli Mlawer
 * email: rrtmgp@aer.com
 *
 * Copyright 2015-2020,  Atmospheric and Environmental Research and
 * Regents of the University of Colorado.  All right reserved.
 *
 * This C++ interface can be downloaded from https://github.com/microhh/rte-rrtmgp-cpp
 *
 * Contact: Chiel van Heerwaarden
 * email: chiel.vanheerwaarden@wur.nl
 *
 * Copyright 2020, Wageningen University & Research.
 *
 * Use and duplication is permitted under the terms of the
 * BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
 *
 */

#include "Rte_lw.h"
#include "Array.h"
#include "Optical_props.h"
#include "Source_functions.h"
#include "Fluxes.h"

#include "rrtmgp_kernels.h"

namespace rrtmgp_kernel_launcher
{
    template<typename TF>
    void apply_BC(
            int ncol, int nlay, int ngpt,
            BOOL_TYPE top_at_1, Array<TF,3>& gpt_flux_dn)
    {
        rrtmgp_kernels::apply_BC_0(
                &ncol, &nlay, &ngpt,
                &top_at_1, gpt_flux_dn.ptr());
    }

    template<typename TF>
    void apply_BC(
            int ncol, int nlay, int ngpt,
            BOOL_TYPE top_at_1, const Array<TF,2>& inc_flux,
            Array<TF,3>& gpt_flux_dn)
    {
        rrtmgp_kernels::apply_BC_gpt(
                &ncol, &nlay, &ngpt,
                &top_at_1, const_cast<TF*>(inc_flux.ptr()), gpt_flux_dn.ptr());
    }

    template<typename TF>
    void lw_solver_noscat_GaussQuad(
            int ncol, int nlay, int ngpt, BOOL_TYPE top_at_1, int n_quad_angs,
            const Array<TF,2>& gauss_Ds_subset,
            const Array<TF,2>& gauss_wts_subset,
            const Array<TF,3>& tau,
            const Array<TF,3>& lay_source,
            const Array<TF,3>& lev_source_inc, const Array<TF,3>& lev_source_dec,
            const Array<TF,2>& sfc_emis_gpt, const Array<TF,2>& sfc_source,
            Array<TF,3>& gpt_flux_up, Array<TF,3>& gpt_flux_dn,
            Array<TF,2>& sfc_source_jac, Array<TF,3>& gpt_flux_up_jac)
{
    rrtmgp_kernels::lw_solver_noscat_GaussQuad(
                &ncol, &nlay, &ngpt, &top_at_1, &n_quad_angs,
                const_cast<TF*>(gauss_Ds_subset.ptr()),
                const_cast<TF*>(gauss_wts_subset.ptr()),
                const_cast<TF*>(tau.ptr()),
                const_cast<TF*>(lay_source.ptr()),
                const_cast<TF*>(lev_source_inc.ptr()),
                const_cast<TF*>(lev_source_dec.ptr()),
                const_cast<TF*>(sfc_emis_gpt.ptr()),
                const_cast<TF*>(sfc_source.ptr()),
                gpt_flux_up.ptr(),
                gpt_flux_dn.ptr(),
                sfc_source_jac.ptr(),
                gpt_flux_up_jac.ptr());
    }
}

template<typename TF>
void Rte_lw<TF>::rte_lw(
        const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
        const BOOL_TYPE top_at_1,
        const Source_func_lw<TF>& sources,
        const Array<TF,2>& sfc_emis,
        const Array<TF,2>& inc_flux,
        Array<TF,3>& gpt_flux_up,
        Array<TF,3>& gpt_flux_dn,
        const int n_gauss_angles)
{
    const int max_gauss_pts = 4;
    const Array<TF,2> gauss_Ds(
            {      1.66,         0.,         0.,         0.,
             1.18350343, 2.81649655,         0.,         0.,
             1.09719858, 1.69338507, 4.70941630,         0.,
             1.06056257, 1.38282560, 2.40148179, 7.15513024},
            { max_gauss_pts, max_gauss_pts });

    const Array<TF,2> gauss_wts(
            {         0.5,           0.,           0.,           0.,
             0.3180413817, 0.1819586183,           0.,           0.,
             0.2009319137, 0.2292411064, 0.0698269799,           0.,
             0.1355069134, 0.2034645680, 0.1298475476, 0.0311809710},
            { max_gauss_pts, max_gauss_pts });

    const int ncol = optical_props->get_ncol();
    const int nlay = optical_props->get_nlay();
    const int ngpt = optical_props->get_ngpt();

    Array<TF,2> sfc_emis_gpt({ncol, ngpt});

    expand_and_transpose(optical_props, sfc_emis, sfc_emis_gpt);

    // Upper boundary condition.
    if (inc_flux.size() == 0)
        rrtmgp_kernel_launcher::apply_BC(ncol, nlay, ngpt, top_at_1, gpt_flux_dn);
    else
        rrtmgp_kernel_launcher::apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, gpt_flux_dn);

    // Run the radiative transfer solver
    const int n_quad_angs = n_gauss_angles;

    Array<TF,2> gauss_Ds_subset = gauss_Ds.subset(
            {{ {1, n_quad_angs}, {n_quad_angs, n_quad_angs} }});
    Array<TF,2> gauss_wts_subset = gauss_wts.subset(
            {{ {1, n_quad_angs}, {n_quad_angs, n_quad_angs} }});

    // For now, just pass the arrays around.
    Array<TF,2> sfc_src_jac(sources.get_sfc_source().get_dims());
    Array<TF,3> gpt_flux_up_jac(gpt_flux_up.get_dims());

    rrtmgp_kernel_launcher::lw_solver_noscat_GaussQuad(
            ncol, nlay, ngpt, top_at_1, n_quad_angs,
            gauss_Ds_subset, gauss_wts_subset,
            optical_props->get_tau(),
            sources.get_lay_source(),
            sources.get_lev_source_inc(), sources.get_lev_source_dec(),
            sfc_emis_gpt, sources.get_sfc_source(),
            gpt_flux_up, gpt_flux_dn,
            sfc_src_jac, gpt_flux_up_jac);

    // CvH: In the fortran code this call is here, I removed it for performance and flexibility.
    // fluxes->reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1);
}

template<typename TF>
void Rte_lw<TF>::expand_and_transpose(
        const std::unique_ptr<Optical_props_arry<TF>>& ops,
        const Array<TF,2> arr_in,
        Array<TF,2>& arr_out)
{
    const int ncol = arr_in.dim(2);
    const int nband = ops->get_nband();
    Array<int,2> limits = ops->get_band_lims_gpoint();

    for (int iband=1; iband<=nband; ++iband)
        for (int icol=1; icol<=ncol; ++icol)
            for (int igpt=limits({1, iband}); igpt<=limits({2, iband}); ++igpt)
                arr_out({icol, igpt}) = arr_in({iband, icol});
}

#ifdef FLOAT_SINGLE_RRTMGP
template class Rte_lw<float>;
#else
template class Rte_lw<double>;
#endif
