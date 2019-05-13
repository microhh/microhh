#ifndef RRTMGP_KERNELS_H
#define RRTMGP_KERNELS_H

#ifdef FLOAT_SINGLE
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
}

// Kernels of gas optics.
namespace rrtmgp_kernels
{
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
                int* tropo,
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
            int* minor_scales_with_density_lower,
            int* minor_scales_with_density_upper,
            int* scale_by_complement_lower,
            int* scale_by_complement_upper,
            int* idx_minor_lower,
            int* idx_minor_upper,
            int* idx_minor_scaling_lower,
            int* idx_minor_scaling_upper,
            int* kminor_start_lower,
            int* kminor_start_upper,
            int* tropo,
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
            FLOAT_TYPE* fmajor, int* jeta, int* tropo, int* jtemp, int* jpress,
            int* gpoint_bands, int* band_lims_gpt, FLOAT_TYPE* pfracin, FLOAT_TYPE* temp_ref_min,
            FLOAT_TYPE* totplnk_delta, FLOAT_TYPE* totplnk, int* gpoint_flavor,
            FLOAT_TYPE* sfc_src, FLOAT_TYPE* lay_src, FLOAT_TYPE* lev_src, FLOAT_TYPE* lev_source_dec);

    extern "C" void compute_tau_rayleigh(
            int* ncol, int* nlay, int* nband, int* ngpt,
            int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
            int* gpoint_flavor,
            int* band_lims_gpt,
            FLOAT_TYPE* krayl,
            int* idx_h2o, FLOAT_TYPE* col_dry, FLOAT_TYPE* col_gas,
            FLOAT_TYPE* fminor, int* eta,
            int* tropo, int* jtemp,
            FLOAT_TYPE* tau_rayleigh);
}

// Kernels of longwave solver.
namespace rrtmgp_kernels
{
    extern "C" void apply_BC_0(
            int* ncol, int* nlay, int* ngpt,
            int* top_at_1, FLOAT_TYPE* gpt_flux_dn);

    extern "C" void lw_solver_noscat_GaussQuad(
            int* ncol, int* nlay, int* ngpt, int* top_at_1, int* n_quad_angs,
            FLOAT_TYPE* gauss_Ds_subset, FLOAT_TYPE* gauss_wts_subset,
            FLOAT_TYPE* tau,
            FLOAT_TYPE* lay_source, FLOAT_TYPE* lev_source_inc, FLOAT_TYPE* lev_source_dec,
            FLOAT_TYPE* sfc_emis_gpt, FLOAT_TYPE* sfc_source,
            FLOAT_TYPE* gpt_flux_up, FLOAT_TYPE* gpt_flux_dn);
}

// Kernels of shortwave solver.
namespace rrtmgp_kernels
{
    extern "C" void apply_BC_0(
            int* ncol, int* nlay, int* ngpt,
            int* top_at_1, FLOAT_TYPE* gpt_flux_dn);

    extern "C" void apply_BC_factor(
            int* ncol, int* nlay, int* ngpt,
            int* top_at_1, FLOAT_TYPE* inc_flux,
            FLOAT_TYPE* factor, FLOAT_TYPE* flux_dn);

    extern "C" void sw_solver_2stream(
            int* ncol, int* nlay, int* ngpt, int* top_at_1,
            FLOAT_TYPE* tau,
            FLOAT_TYPE* ssa,
            FLOAT_TYPE* g,
            FLOAT_TYPE* mu0,
            FLOAT_TYPE* sfc_alb_dir_gpt, FLOAT_TYPE* sfc_alb_dif_gpt,
            FLOAT_TYPE* gpt_flux_up, FLOAT_TYPE* gpt_flux_dn, FLOAT_TYPE* gpt_flux_dir);
}

#endif
