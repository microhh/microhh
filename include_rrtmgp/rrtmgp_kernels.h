/*
 * This file is part of a C++ interface to the Radiative Transfer for Energetics (RTE)
 * and Rapid Radiative Transfer Model for GCM applications Parallel (RRTMGP).
 *
 * The original code is found at https://github.com/earth-system-radiation/rte-rrtmgp.
 *
 * Contacts: Robert Pincus and Eli Mlawer
 * email: rrtmgp@aer.com
 *
 * Copyright 2015-2020,  Atmospheric and Environmental Research and
 * Regents of the University of Colorado.  All right reserved.
 *
 * This C++ interface can be downloaded from https://github.com/earth-system-radiation/rte-rrtmgp-cpp
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

#ifndef RRTMGP_KERNELS_H
#define RRTMGP_KERNELS_H

#ifdef FLOAT_SINGLE_RRTMGP
#define FLOAT_TYPE float
#else
#define FLOAT_TYPE double
#endif

// Kernels of fluxes.
namespace rrtmgp_kernels
{
    extern "C" void sum_broadband(
            int* ncol, int* nlev, int* ngpt,
            FLOAT_TYPE* spectral_flux, FLOAT_TYPE* broadband_flux);

    extern "C" void net_broadband_precalc(
            int* ncol, int* nlev,
            FLOAT_TYPE* broadband_flux_dn, FLOAT_TYPE* broadband_flux_up,
            FLOAT_TYPE* broadband_flux_net);

    extern "C" void sum_byband(
            int* ncol, int* nlev, int* ngpt, int* nbnd,
            int* band_lims,
            FLOAT_TYPE* spectral_flux,
            FLOAT_TYPE* byband_flux);

    extern "C" void net_byband_precalc(
            int* ncol, int* nlev, int* nbnd,
            FLOAT_TYPE* byband_flux_dn, FLOAT_TYPE* byband_flux_up,
            FLOAT_TYPE* byband_flux_net);

    extern "C" void zero_array_3D(
            int* ni, int* nj, int* nk, FLOAT_TYPE* array);

    extern "C" void zero_array_4D(
             int* ni, int* nj, int* nk, int* nl, FLOAT_TYPE* array);

    extern "C" void interpolation(
                int* ncol, int* nlay,
                int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
                int* flavor,
                FLOAT_TYPE* press_ref_log,
                FLOAT_TYPE* temp_ref,
                FLOAT_TYPE* press_ref_log_delta,
                FLOAT_TYPE* temp_ref_min,
                FLOAT_TYPE* temp_ref_delta,
                FLOAT_TYPE* press_ref_trop_log,
                FLOAT_TYPE* vmr_ref,
                FLOAT_TYPE* play,
                FLOAT_TYPE* tlay,
                FLOAT_TYPE* col_gas,
                int* jtemp,
                FLOAT_TYPE* fmajor, FLOAT_TYPE* fminor,
                FLOAT_TYPE* col_mix,
                BOOL_TYPE* tropo,
                int* jeta,
                int* jpress);

    extern "C" void compute_tau_absorption(
            int* ncol, int* nlay, int* nband, int* ngpt,
            int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
            int* nminorlower, int* nminorklower,
            int* nminorupper, int* nminorkupper,
            int* idx_h2o,
            int* gpoint_flavor,
            int* band_lims_gpt,
            FLOAT_TYPE* kmajor,
            FLOAT_TYPE* kminor_lower,
            FLOAT_TYPE* kminor_upper,
            int* minor_limits_gpt_lower,
            int* minor_limits_gpt_upper,
            BOOL_TYPE* minor_scales_with_density_lower,
            BOOL_TYPE* minor_scales_with_density_upper,
            BOOL_TYPE* scale_by_complement_lower,
            BOOL_TYPE* scale_by_complement_upper,
            int* idx_minor_lower,
            int* idx_minor_upper,
            int* idx_minor_scaling_lower,
            int* idx_minor_scaling_upper,
            int* kminor_start_lower,
            int* kminor_start_upper,
            BOOL_TYPE* tropo,
            FLOAT_TYPE* col_mix, FLOAT_TYPE* fmajor, FLOAT_TYPE* fminor,
            FLOAT_TYPE* play, FLOAT_TYPE* tlay, FLOAT_TYPE* col_gas,
            int* jeta, int* jtemp, int* jpress,
            FLOAT_TYPE* tau);

    extern "C" void reorder_123x321_kernel(
            int* dim1, int* dim2, int* dim3,
            FLOAT_TYPE* array, FLOAT_TYPE* array_out);

    extern "C" void combine_and_reorder_2str(
            int* ncol, int* nlay, int* ngpt,
            FLOAT_TYPE* tau_local, FLOAT_TYPE* tau_rayleigh,
            FLOAT_TYPE* tau, FLOAT_TYPE* ssa, FLOAT_TYPE* g);

    extern "C" void compute_Planck_source(
            int* ncol, int* nlay, int* nbnd, int* ngpt,
            int* nflav, int* neta, int* npres, int* ntemp, int* nPlanckTemp,
            FLOAT_TYPE* tlay, FLOAT_TYPE* tlev, FLOAT_TYPE* tsfc, int* sfc_lay,
            FLOAT_TYPE* fmajor, int* jeta, BOOL_TYPE* tropo, int* jtemp, int* jpress,
            int* gpoint_bands, int* band_lims_gpt, FLOAT_TYPE* pfracin, FLOAT_TYPE* temp_ref_min,
            FLOAT_TYPE* totplnk_delta, FLOAT_TYPE* totplnk, int* gpoint_flavor,
            FLOAT_TYPE* sfc_src, FLOAT_TYPE* lay_src, FLOAT_TYPE* lev_src, FLOAT_TYPE* lev_source_dec,
            FLOAT_TYPE* sfc_src_jac);

    extern "C" void compute_tau_rayleigh(
            int* ncol, int* nlay, int* nband, int* ngpt,
            int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
            int* gpoint_flavor,
            int* band_lims_gpt,
            FLOAT_TYPE* krayl,
            int* idx_h2o, FLOAT_TYPE* col_dry, FLOAT_TYPE* col_gas,
            FLOAT_TYPE* fminor, int* eta,
            BOOL_TYPE* tropo, int* jtemp,
            FLOAT_TYPE* tau_rayleigh);

    extern "C" void apply_BC_0(
            int* ncol, int* nlay, int* ngpt,
            BOOL_TYPE* top_at_1, FLOAT_TYPE* gpt_flux_dn);

    extern "C" void apply_BC_gpt(
            int* ncol, int* nlay, int* ngpt,
            BOOL_TYPE* top_at_1, FLOAT_TYPE* inc_flux, FLOAT_TYPE* gpt_flux_dn);

    extern "C" void lw_solver_noscat_GaussQuad(
            int* ncol, int* nlay, int* ngpt, BOOL_TYPE* top_at_1, int* n_quad_angs,
            FLOAT_TYPE* gauss_Ds_subset, FLOAT_TYPE* gauss_wts_subset,
            FLOAT_TYPE* tau,
            FLOAT_TYPE* lay_source, FLOAT_TYPE* lev_source_inc, FLOAT_TYPE* lev_source_dec,
            FLOAT_TYPE* sfc_emis_gpt, FLOAT_TYPE* sfc_source,
            FLOAT_TYPE* gpt_flux_up, FLOAT_TYPE* gpt_flux_dn,
            FLOAT_TYPE* sfc_source_jac, FLOAT_TYPE* gpt_flux_up_jac);

    extern "C" void apply_BC_factor(
            int* ncol, int* nlay, int* ngpt,
            BOOL_TYPE* top_at_1, FLOAT_TYPE* inc_flux,
            FLOAT_TYPE* factor, FLOAT_TYPE* flux_dn);

    extern "C" void sw_solver_2stream(
            int* ncol, int* nlay, int* ngpt, BOOL_TYPE* top_at_1,
            FLOAT_TYPE* tau,
            FLOAT_TYPE* ssa,
            FLOAT_TYPE* g,
            FLOAT_TYPE* mu0,
            FLOAT_TYPE* sfc_alb_dir_gpt, FLOAT_TYPE* sfc_alb_dif_gpt,
            FLOAT_TYPE* gpt_flux_up, FLOAT_TYPE* gpt_flux_dn, FLOAT_TYPE* gpt_flux_dir);

    extern "C" void increment_2stream_by_2stream(
            int* ncol, int* nlev, int* ngpt,
            FLOAT_TYPE* tau_inout, FLOAT_TYPE* ssa_inout, FLOAT_TYPE* g_inout,
            FLOAT_TYPE* tau_in, FLOAT_TYPE* ssa_in, FLOAT_TYPE* g_in);

    extern "C" void increment_1scalar_by_1scalar(
            int* ncol, int* nlev, int* ngpt,
            FLOAT_TYPE* tau_inout, FLOAT_TYPE* tau_in);

    extern "C" void inc_2stream_by_2stream_bybnd(
            int* ncol, int* nlev, int* ngpt,
            FLOAT_TYPE* tau_inout, FLOAT_TYPE* ssa_inout, FLOAT_TYPE* g_inout,
            FLOAT_TYPE* tau_in, FLOAT_TYPE* ssa_in, FLOAT_TYPE* g_in,
            int* nbnd, int* band_lims_gpoint);

    extern "C" void inc_1scalar_by_1scalar_bybnd(
            int* ncol, int* nlev, int* ngpt,
            FLOAT_TYPE* tau_inout, FLOAT_TYPE* tau_in,
            int* nbnd, int* band_lims_gpoint);

    extern "C" void delta_scale_2str_k(
            int* ncol, int* nlev, int* ngpt,
            FLOAT_TYPE* tau_inout, FLOAT_TYPE* ssa_inout, FLOAT_TYPE* g_inout);
}
#endif
