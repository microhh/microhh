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

#include "Rte_sw.h"
#include "Array.h"
#include "Optical_props.h"
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
            int ncol, int nlay, int ngpt, BOOL_TYPE top_at_1,
            const Array<TF,2>& inc_flux, Array<TF,3>& gpt_flux_dn)
    {
        rrtmgp_kernels::apply_BC_gpt(
                &ncol, &nlay, &ngpt, &top_at_1,
                const_cast<TF*>(inc_flux.ptr()), gpt_flux_dn.ptr());
    }

    template<typename TF>
    void apply_BC(
            int ncol, int nlay, int ngpt, BOOL_TYPE top_at_1,
            const Array<TF,2>& inc_flux,
            const Array<TF,1>& factor,
            Array<TF,3>& gpt_flux)
    {
        rrtmgp_kernels::apply_BC_factor(
                &ncol, &nlay, &ngpt,
                &top_at_1,
                const_cast<TF*>(inc_flux.ptr()),
                const_cast<TF*>(factor.ptr()),
                gpt_flux.ptr());
    }

    template<typename TF>
    void sw_solver_2stream(
            int ncol, int nlay, int ngpt, BOOL_TYPE top_at_1,
            const Array<TF,3>& tau,
            const Array<TF,3>& ssa,
            const Array<TF,3>& g,
            const Array<TF,1>& mu0,
            const Array<TF,2>& sfc_alb_dir_gpt, const Array<TF,2>& sfc_alb_dif_gpt,
            Array<TF,3>& gpt_flux_up, Array<TF,3>& gpt_flux_dn, Array<TF,3>& gpt_flux_dir)
    {
        rrtmgp_kernels::sw_solver_2stream(
                &ncol, &nlay, &ngpt, &top_at_1,
                const_cast<TF*>(tau.ptr()),
                const_cast<TF*>(ssa.ptr()),
                const_cast<TF*>(g  .ptr()),
                const_cast<TF*>(mu0.ptr()),
                const_cast<TF*>(sfc_alb_dir_gpt.ptr()),
                const_cast<TF*>(sfc_alb_dif_gpt.ptr()),
                gpt_flux_up.ptr(), gpt_flux_dn.ptr(), gpt_flux_dir.ptr());
    }
}

template<typename TF>
void Rte_sw<TF>::rte_sw(
        const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
        const BOOL_TYPE top_at_1,
        const Array<TF,1>& mu0,
        const Array<TF,2>& inc_flux_dir,
        const Array<TF,2>& sfc_alb_dir,
        const Array<TF,2>& sfc_alb_dif,
        const Array<TF,2>& inc_flux_dif,
        Array<TF,3>& gpt_flux_up,
        Array<TF,3>& gpt_flux_dn,
        Array<TF,3>& gpt_flux_dir)
{
    const int ncol = optical_props->get_ncol();
    const int nlay = optical_props->get_nlay();
    const int ngpt = optical_props->get_ngpt();

    Array<TF,2> sfc_alb_dir_gpt({ncol, ngpt});
    Array<TF,2> sfc_alb_dif_gpt({ncol, ngpt});

    expand_and_transpose(optical_props, sfc_alb_dir, sfc_alb_dir_gpt);
    expand_and_transpose(optical_props, sfc_alb_dif, sfc_alb_dif_gpt);

    // Upper boundary condition. At this stage, flux_dn contains the diffuse radiation only.
    rrtmgp_kernel_launcher::apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dir, mu0, gpt_flux_dir);
    if (inc_flux_dif.size() == 0)
        rrtmgp_kernel_launcher::apply_BC(ncol, nlay, ngpt, top_at_1, gpt_flux_dn);
    else
        rrtmgp_kernel_launcher::apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif, gpt_flux_dn);

    // Run the radiative transfer solver
    // CvH: only two-stream solutions, I skipped the sw_solver_noscat
    rrtmgp_kernel_launcher::sw_solver_2stream(
            ncol, nlay, ngpt, top_at_1,
            optical_props->get_tau(),
            optical_props->get_ssa(),
            optical_props->get_g  (),
            mu0,
            sfc_alb_dir_gpt, sfc_alb_dif_gpt,
            gpt_flux_up, gpt_flux_dn, gpt_flux_dir);

    // CvH: The original fortran code had a call to the reduce here.
    // fluxes->reduce(gpt_flux_up, gpt_flux_dn, gpt_flux_dir, optical_props, top_at_1);
}

template<typename TF>
void Rte_sw<TF>::expand_and_transpose(
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
template class Rte_sw<float>;
#else
template class Rte_sw<double>;
#endif
