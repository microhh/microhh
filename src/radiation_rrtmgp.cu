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

#include "radiation_rrtmgp.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "thermo.h"
#include "stats.h"
#include "constants.h"

#include "Array.h"


namespace
{
    #ifdef USECUDA

    __global__
    void calc_tendency(
            Float* __restrict__ thlt_rad,  const Float* __restrict__ flux_up, 
            const Float* __restrict flux_dn, const Float* __restrict__ rho, 
            const Float* __restrict__ exner, const Float* __restrict__ dz,
            const int istart, const int jstart, const int kstart,
            const int iend,   const int jend,   const int kend,
            const int igc, const int jgc, const int kgc,
            const int jj, const int kk,
            const int jj_nogc, const int kk_nogc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z*blockDim.z + threadIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const Float fac = Float(1.) / (rho[k] * Constants::cp<Float> * exner[k] * dz[k]);

            const int ijk = i + j*jj + k*kk;
            const int ijk_nogc = (i-igc) + (j-jgc)*jj_nogc + (k-kgc)*kk_nogc;

            thlt_rad[ijk] = fac * (flux_up[ijk_nogc + kk_nogc] - flux_up[ijk_nogc] -
                                   flux_dn[ijk_nogc + kk_nogc] + flux_dn[ijk_nogc] ); 
        }
    }
            
    __global__
    void add_tendency(
            Float* __restrict__ thlt,  const Float* __restrict__ thlt_rad,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z*blockDim.z + threadIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            thlt[ijk] = thlt_rad[ijk];
        }
    }
    #endif
}


#ifdef USECUDA
template<typename TF>
void Radiation_rrtmgp<TF>::exec_longwave(
        Thermo<TF>& thermo, Timeloop<TF>& timeloop, Stats<TF>& stats,
        Array_gpu<Float,2>& flux_up, Array_gpu<Float,2>& flux_dn, Array_gpu<Float,2>& flux_net,
        const Array_gpu<Float,2>& t_lay, const Array_gpu<Float,2>& t_lev, const Array_gpu<Float,1>& t_sfc,
        const Array_gpu<Float,2>& h2o, const Array_gpu<Float,2>& clwp, const Array_gpu<Float,2>& ciwp,
        const bool compute_clouds)
{
    constexpr int n_col_block = 1024;

    auto& gd = grid.get_grid_data();

    const int n_col = gd.imax*gd.jmax;
    const int n_lay = gd.ktot;
    const int n_lev = gd.ktot+1;

    const int n_blocks = n_col / n_col_block;
    const int n_col_block_left = n_col % n_col_block;

    const int n_gpt = this->kdist_lw_gpu->get_ngpt();
    const int n_bnd = this->kdist_lw_gpu->get_nband();

    const Bool top_at_1 = 0;

    // Define the pointers for the subsetting.
    std::unique_ptr<Optical_props_arry_gpu> optical_props_subset =
            std::make_unique<Optical_props_1scl_gpu>(n_col_block, n_lay, *kdist_lw_gpu);
    std::unique_ptr<Source_func_lw_gpu> sources_subset =
            std::make_unique<Source_func_lw_gpu>(n_col_block, n_lay, *kdist_lw_gpu);
    std::unique_ptr<Optical_props_1scl_gpu> cloud_optical_props_subset =
            std::make_unique<Optical_props_1scl_gpu>(n_col_block, n_lay, *cloud_lw_gpu);

    std::unique_ptr<Optical_props_arry_gpu> optical_props_left =
            std::make_unique<Optical_props_1scl_gpu>(n_col_block_left, n_lay, *kdist_lw_gpu);
    std::unique_ptr<Source_func_lw_gpu> sources_left =
            std::make_unique<Source_func_lw_gpu>(n_col_block_left, n_lay, *kdist_lw_gpu);
    std::unique_ptr<Optical_props_1scl_gpu> cloud_optical_props_left =
            std::make_unique<Optical_props_1scl_gpu>(n_col_block_left, n_lay, *cloud_lw_gpu);

    /*
    // Lambda function for solving optical properties subset.
    auto call_kernels = [&](
            const int col_s_in, const int col_e_in,
            std::unique_ptr<Optical_props_arry_gpu>& optical_props_subset_in,
            std::unique_ptr<Optical_props_1scl_gpu>& cloud_optical_props_subset_in,
            Source_func_lw_gpu& sources_subset_in,
            const Array_gpu<Float,2>& emis_sfc_subset_in,
            Fluxes_broadband_gpu& fluxes,
            Fluxes_broadband_gpu& bnd_fluxes)
    {
        const int n_col_in = col_e_in - col_s_in + 1;
        Gas_concs_gpu gas_concs_subset(gas_concs, col_s_in, n_col_in);

        auto p_lev_subset = p_lev.subset({{ {col_s_in, col_e_in}, {1, n_lev} }});

        Array_gpu<Float,2> col_dry_subset({n_col_in, n_lay});
        if (col_dry.size() == 0)
            Gas_optics_rrtmgp_gpu::get_col_dry(col_dry_subset, gas_concs_subset.get_vmr("h2o"), p_lev_subset);
        else
            col_dry_subset = col_dry.subset({{ {col_s_in, col_e_in}, {1, n_lay} }});

        kdist_gpu->gas_optics(
                p_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                p_lev_subset,
                t_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                t_sfc.subset({{ {col_s_in, col_e_in} }}),
                gas_concs_subset,
                optical_props_subset_in,
                sources_subset_in,
                col_dry_subset,
                t_lev.subset({{ {col_s_in, col_e_in}, {1, n_lev} }}) );

        if (switch_cloud_optics)
        {
            cloud_optics_gpu->cloud_optics(
                    lwp.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                    iwp.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                    rel.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                    rei.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                    *cloud_optical_props_subset_in);

            // cloud->delta_scale();

            // Add the cloud optical props to the gas optical properties.
            add_to(
                    dynamic_cast<Optical_props_1scl_gpu&>(*optical_props_subset_in),
                    dynamic_cast<Optical_props_1scl_gpu&>(*cloud_optical_props_subset_in));
        }

        // Store the optical properties, if desired.
        if (switch_output_optical)
        {
            subset_kernel_launcher_cuda::get_from_subset(
                    n_col, n_lay, n_gpt, n_col_in, col_s_in, tau.ptr(), lay_source.ptr(), lev_source_inc.ptr(), lev_source_dec.ptr(),
                    optical_props_subset_in->get_tau().ptr(), sources_subset_in.get_lay_source().ptr(),
                    sources_subset_in.get_lev_source_inc().ptr(), sources_subset_in.get_lev_source_dec().ptr());

            subset_kernel_launcher_cuda::get_from_subset(
                    n_col, n_gpt, n_col_in, col_s_in, sfc_source.ptr(), sources_subset_in.get_sfc_source().ptr());
        }

        if (!switch_fluxes)
            return;

        Array_gpu<Float,3> gpt_flux_up({n_col_in, n_lev, n_gpt});
        Array_gpu<Float,3> gpt_flux_dn({n_col_in, n_lev, n_gpt});

        constexpr int n_ang = 1;

        rte_lw.rte_lw(
                optical_props_subset_in,
                top_at_1,
                sources_subset_in,
                emis_sfc_subset_in,
                Array_gpu<Float,2>(), // Add an empty array, no inc_flux.
                gpt_flux_up, gpt_flux_dn,
                n_ang);

        fluxes.reduce(gpt_flux_up, gpt_flux_dn, optical_props_subset_in, top_at_1);

        // Copy the data to the output.
        subset_kernel_launcher_cuda::get_from_subset(
                n_col, n_lev, n_col_in, col_s_in, lw_flux_up.ptr(), lw_flux_dn.ptr(), lw_flux_net.ptr(),
                fluxes.get_flux_up().ptr(), fluxes.get_flux_dn().ptr(), fluxes.get_flux_net().ptr());

        if (switch_output_bnd_fluxes)
        {
            bnd_fluxes.reduce(gpt_flux_up, gpt_flux_dn, optical_props_subset_in, top_at_1);

            subset_kernel_launcher_cuda::get_from_subset(
                    n_col, n_lev, n_bnd, n_col_in, col_s_in, lw_bnd_flux_up.ptr(), lw_bnd_flux_dn.ptr(), lw_bnd_flux_net.ptr(),
                    bnd_fluxes.get_bnd_flux_up().ptr(), bnd_fluxes.get_bnd_flux_dn().ptr(), bnd_fluxes.get_bnd_flux_net().ptr());
        }
    };

    for (int b=1; b<=n_blocks; ++b)
    {
        const int col_s = (b-1) * n_col_block + 1;
        const int col_e =  b    * n_col_block;

        Array_gpu<Float,2> emis_sfc_subset = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});

        std::unique_ptr<Fluxes_broadband_gpu> fluxes_subset =
                std::make_unique<Fluxes_broadband_gpu>(n_col_block, n_lev);
        std::unique_ptr<Fluxes_broadband_gpu> bnd_fluxes_subset =
                std::make_unique<Fluxes_byband_gpu>(n_col_block, n_lev, n_bnd);

        call_kernels(
                col_s, col_e,
                optical_props_subset,
                cloud_optical_props_subset,
                *sources_subset,
                emis_sfc_subset,
                *fluxes_subset,
                *bnd_fluxes_subset);
    }

    if (n_col_block_residual > 0)
    {
        const int col_s = n_col - n_col_block_residual + 1;
        const int col_e = n_col;

        Array_gpu<Float,2> emis_sfc_residual = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});
        std::unique_ptr<Fluxes_broadband_gpu> fluxes_residual =
                std::make_unique<Fluxes_broadband_gpu>(n_col_block_residual, n_lev);
        std::unique_ptr<Fluxes_broadband_gpu> bnd_fluxes_residual =
                std::make_unique<Fluxes_byband_gpu>(n_col_block_residual, n_lev, n_bnd);

        call_kernels(
                col_s, col_e,
                optical_props_residual,
                cloud_optical_props_residual,
                *sources_residual,
                emis_sfc_residual,
                *fluxes_residual,
                *bnd_fluxes_residual);
    }
    */
}
#endif


#ifdef USECUDA
template <typename TF>
void Radiation_rrtmgp<TF>::exec(Thermo<TF>& thermo, double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU_3d (gridi, gridj, gd.kmax+1);
    dim3 blockGPU_3d(blocki, blockj, 1);
    
    const bool do_radiation = ((timeloop.get_itime() % idt_rad == 0) && !timeloop.in_substep()) ;

    if (do_radiation)
    {
        // Set the tendency to zero.
        std::fill(fields.sd.at("thlt_rad")->fld.begin(), fields.sd.at("thlt_rad")->fld.end(), Float(0.));

        auto t_lay = fields.get_tmp_g();
        auto t_lev = fields.get_tmp_g();
        auto h2o   = fields.get_tmp_g(); // This is the volume mixing ratio, not the specific humidity of vapor.
        auto clwp  = fields.get_tmp_g();
        auto ciwp  = fields.get_tmp_g();

        // Set the input to the radiation on a 3D grid without ghost cells.
        thermo.get_radiation_fields_g(*t_lay, *t_lev, *h2o, *clwp, *ciwp);

        const int nmaxh = gd.imax*gd.jmax*(gd.ktot+1);
        const int ijmax = gd.imax*gd.jmax;

        // Create views on existing variables.
        Array_gpu<Float,2> t_lay_a(t_lay->fld_g, {gd.imax*gd.jmax, gd.ktot});
        Array_gpu<Float,2> t_lev_a(t_lev->fld_g, {gd.imax*gd.jmax, gd.ktot+1});
        Array_gpu<Float,1> t_sfc_a(t_lev->fld_bot_g, {gd.imax*gd.jmax});
        Array_gpu<Float,2> h2o_a(h2o->fld_g, {gd.imax*gd.jmax, gd.ktot});
        Array_gpu<Float,2> clwp_a(clwp->fld_g, {gd.imax*gd.jmax, gd.ktot});
        Array_gpu<Float,2> ciwp_a(ciwp->fld_g, {gd.imax*gd.jmax, gd.ktot});

        // Flux fields.
        Array_gpu<Float,2> flux_up ({gd.imax*gd.jmax, gd.ktot+1});
        Array_gpu<Float,2> flux_dn ({gd.imax*gd.jmax, gd.ktot+1});
        Array_gpu<Float,2> flux_net({gd.imax*gd.jmax, gd.ktot+1});

        const bool compute_clouds = true;
 
        try
        {
            if (sw_longwave)
            {
                /*
                exec_longwave(
                        thermo, timeloop, stats,
                        flux_up, flux_dn, flux_net,
                        t_lay_a, t_lev_a, t_sfc_a, h2o_a, clwp_a, ciwp_a,
                        compute_clouds); */
                calc_tendency<<<gridGPU_3d, blockGPU_3d>>>(
                        fields.sd.at("thlt_rad")->fld.data(),
                        flux_up.ptr(), flux_dn.ptr(),
                        fields.rhoref.data(), thermo.get_basestate_vector("exner").data(),
                        gd.dz.data(),
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.igc, gd.jgc, gd.kgc,
                        gd.icells, gd.ijcells,
                        gd.imax, gd.imax*gd.jmax);
                /*
                store_surface_fluxes(
                        lw_flux_up_sfc.data(), lw_flux_dn_sfc.data(),
                        flux_up.ptr(), flux_dn.ptr(),
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.igc, gd.jgc,
                        gd.icells, gd.ijcells,
                        gd.imax);
                        */
            }

            /*
            if (sw_shortwave)
            {
                if (!sw_fixed_sza)
                {
                    // Update the solar zenith angle, and calculate new shortwave reference column
                    const int day_of_year = int(timeloop.calc_day_of_year());
                    const int year = timeloop.get_year();
                    const Float seconds_after_midnight = Float(timeloop.calc_hour_of_day()*3600);
                    this->mu0 = calc_cos_zenith_angle(lat, lon, day_of_year, seconds_after_midnight, year);

                    // Calculate correction factor for impact Sun's distance on the solar "constant"
                    const Float frac_day_of_year = Float(day_of_year) + seconds_after_midnight / Float(86400);
                    this->tsi_scaling = calc_sun_distance_factor(frac_day_of_year);

                    if (is_day(this->mu0))
                    {
                        const int n_bnd = kdist_sw->get_nband();

                        // Set the solar zenith angle and albedo.
                        Array<Float,2> sfc_alb_dir({n_bnd, n_col});
                        Array<Float,2> sfc_alb_dif({n_bnd, n_col});

                        for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
                        {
                            sfc_alb_dir({ibnd, 1}) = this->sfc_alb_dir;
                            sfc_alb_dif({ibnd, 1}) = this->sfc_alb_dif;
                        }

                        Array<Float,1> mu0({n_col});
                        mu0({1}) = this->mu0;

                        solve_shortwave_column(
                                optical_props_sw,
                                sw_flux_up_col, sw_flux_dn_col, sw_flux_dn_dir_col, sw_flux_net_col,
                                sw_flux_dn_dir_inc, sw_flux_dn_dif_inc, thermo.get_basestate_vector("ph")[gd.kend],
                                gas_concs_col,
                                *kdist_sw,
                                col_dry,
                                p_lay_col, p_lev_col,
                                t_lay_col, t_lev_col,
                                mu0,
                                sfc_alb_dir, sfc_alb_dif,
                                tsi_scaling,
                                n_lay_col);
                    }
                }

                if (is_day(this->mu0))
                {
                    Array<Float,2> flux_dn_dir({gd.imax*gd.jmax, gd.ktot+1});

                    exec_shortwave(
                            thermo, timeloop, stats,
                            flux_up, flux_dn, flux_dn_dir, flux_net,
                            t_lay_a, t_lev_a, h2o_a, clwp_a, ciwp_a,
                            compute_clouds);

                    calc_tendency(
                            fields.sd.at("thlt_rad")->fld.data(),
                            flux_up.ptr(), flux_dn.ptr(),
                            fields.rhoref.data(), thermo.get_basestate_vector("exner").data(),
                            gd.dz.data(),
                            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                            gd.igc, gd.jgc, gd.kgc,
                            gd.icells, gd.ijcells,
                            gd.imax, gd.imax*gd.jmax);

                    store_surface_fluxes(
                            sw_flux_up_sfc.data(), sw_flux_dn_sfc.data(),
                            flux_up.ptr(), flux_dn.ptr(),
                            gd.istart, gd.iend,
                            gd.jstart, gd.jend,
                            gd.igc, gd.jgc,
                            gd.icells, gd.ijcells,
                            gd.imax);
                }
                else
                {
                    // Set the surface fluxes to zero, for (e.g.) the land-surface model.
                    std::fill(sw_flux_up_sfc.begin(), sw_flux_up_sfc.end(), Float(0));
                    std::fill(sw_flux_dn_sfc.begin(), sw_flux_dn_sfc.end(), Float(0));
                }
            }
            */
        } // End try block.
        catch (std::exception& e)
        {
            #ifdef USEMPI
            std::cout << "SINGLE PROCESS EXCEPTION: " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            #else
            throw;
            #endif
        }

        fields.release_tmp_g(t_lay);
        fields.release_tmp_g(t_lev);
        fields.release_tmp_g(h2o);
        fields.release_tmp_g(clwp);
        fields.release_tmp_g(ciwp);
    }

    // Always add the tendency.
    add_tendency<<<gridGPU_3d, blockGPU_3d>>>(
            fields.st.at("thl")->fld_g,
            fields.sd.at("thlt_rad")->fld_g,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_tend(*fields.st.at("thl"), tend_name);
}


template <typename TF>
std::vector<TF>& Radiation_rrtmgp<TF>::get_surface_radiation(const std::string& name)
{
    throw std::runtime_error("Radiation_rrtmgp is not implemented yet on the GPU");
}


template <typename TF>
void Radiation_rrtmgp<TF>::prepare_device()
{
}


template <typename TF>
void Radiation_rrtmgp<TF>::clear_device()
{
}
#endif


#ifdef FLOAT_SINGLE
template class Radiation_rrtmgp<float>;
#else
template class Radiation_rrtmgp<double>;
#endif
