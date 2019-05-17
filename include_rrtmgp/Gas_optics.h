#ifndef GAS_OPTICS_H
#define GAS_OPTICS_H

#include <string>

#include "Array.h"

// Forward declarations.
template<typename TF> class Optical_props;
template<typename TF> class Optical_props_arry;
template<typename TF> class Gas_concs;
template<typename TF> class Source_func_lw;

template<typename TF>
class Gas_optics : public Optical_props<TF>
{
    public:
        // Constructor for longwave variant.
        Gas_optics(
                const Gas_concs<TF>& available_gases,
                const Array<std::string,1>& gas_names,
                const Array<int,3>& key_species,
                const Array<int,2>& band2gpt,
                const Array<TF,2>& band_lims_wavenum,
                const Array<TF,1>& press_ref,
                const TF press_ref_trop,
                const Array<TF,1>& temp_ref,
                const TF temp_ref_p,
                const TF temp_ref_t,
                const Array<TF,3>& vmr_ref,
                const Array<TF,4>& kmajor,
                const Array<TF,3>& kminor_lower,
                const Array<TF,3>& kminor_upper,
                const Array<std::string,1>& gas_minor,
                const Array<std::string,1>& identifier_minor,
                const Array<std::string,1>& minor_gases_lower,
                const Array<std::string,1>& minor_gases_upper,
                const Array<int,2>& minor_limits_gpt_lower,
                const Array<int,2>& minor_limits_gpt_upper,
                const Array<int,1>& minor_scales_with_density_lower,
                const Array<int,1>& minor_scales_with_density_upper,
                const Array<std::string,1>& scaling_gas_lower,
                const Array<std::string,1>& scaling_gas_upper,
                const Array<int,1>& scale_by_complement_lower,
                const Array<int,1>& scale_by_complement_upper,
                const Array<int,1>& kminor_start_lower,
                const Array<int,1>& kminor_start_upper,
                const Array<TF,2>& totplnk,
                const Array<TF,4>& planck_frac,
                const Array<TF,3>& rayl_lower,
                const Array<TF,3>& rayl_upper);

        // Constructor for longwave variant.
        Gas_optics(
                const Gas_concs<TF>& available_gases,
                const Array<std::string,1>& gas_names,
                const Array<int,3>& key_species,
                const Array<int,2>& band2gpt,
                const Array<TF,2>& band_lims_wavenum,
                const Array<TF,1>& press_ref,
                const TF press_ref_trop,
                const Array<TF,1>& temp_ref,
                const TF temp_ref_p,
                const TF temp_ref_t,
                const Array<TF,3>& vmr_ref,
                const Array<TF,4>& kmajor,
                const Array<TF,3>& kminor_lower,
                const Array<TF,3>& kminor_upper,
                const Array<std::string,1>& gas_minor,
                const Array<std::string,1>& identifier_minor,
                const Array<std::string,1>& minor_gases_lower,
                const Array<std::string,1>& minor_gases_upper,
                const Array<int,2>& minor_limits_gpt_lower,
                const Array<int,2>& minor_limits_gpt_upper,
                const Array<int,1>& minor_scales_with_density_lower,
                const Array<int,1>& minor_scales_with_density_upper,
                const Array<std::string,1>& scaling_gas_lower,
                const Array<std::string,1>& scaling_gas_upper,
                const Array<int,1>& scale_by_complement_lower,
                const Array<int,1>& scale_by_complement_upper,
                const Array<int,1>& kminor_start_lower,
                const Array<int,1>& kminor_start_upper,
                const Array<TF,1>& solar_src,
                const Array<TF,3>& rayl_lower,
                const Array<TF,3>& rayl_upper);

        void get_col_dry(
                Array<TF,2>& col_dry,
                const Array<TF,2>& vmr_h2o,
                const Array<TF,2>& plev);

        bool source_is_internal() const { return (totplnk.size() > 0) && (planck_frac.size() > 0); }
        TF get_press_ref_min() const { return press_ref_min; }

        int get_nflav() const { return flavor.dim(2); }
        int get_neta() const { return kmajor.dim(2); }
        int get_npres() const { return kmajor.dim(3)-1; }
        int get_ntemp() const { return kmajor.dim(4); }
        int get_nPlanckTemp() const { return totplnk.dim(1); }

        // Longwave variant.
        void gas_optics(
                const Array<TF,2>& play,
                const Array<TF,2>& plev,
                const Array<TF,2>& tlay,
                const Array<TF,1>& tsfc,
                const Gas_concs<TF>& gas_desc,
                std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                Source_func_lw<TF>& sources,
                const Array<TF,2>& col_dry,
                const Array<TF,2>& tlev);

        // Shortwave variant.
        void gas_optics(
                const Array<TF,2>& play,
                const Array<TF,2>& plev,
                const Array<TF,2>& tlay,
                const Gas_concs<TF>& gas_desc,
                std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                Array<TF,2>& toa_src,
                const Array<TF,2>& col_dry);

    void combine_and_reorder(
            const Array<TF,3>& tau,
            const Array<TF,3>& tau_rayleigh,
            const bool has_rayleigh,
            std::unique_ptr<Optical_props_arry<TF>>& optical_props);

    void source(
            const int ncol, const int nlay, const int nband, const int ngpt,
            const Array<TF,2>& play, const Array<TF,2>& plev,
            const Array<TF,2>& tlay, const Array<TF,1>& tsfc,
            const Array<int,2>& jtemp, const Array<int,2>& jpress,
            const Array<int,4>& jeta, const Array<int,2>& tropo,
            const Array<TF,6>& fmajor,
            Source_func_lw<TF>& sources,
            const Array<TF,2>& tlev);

    private:
        Array<TF,2> totplnk;
        Array<TF,4> planck_frac;
        TF totplnk_delta;
        TF temp_ref_min, temp_ref_max;
        TF press_ref_min, press_ref_max;
        TF press_ref_trop_log;

        TF press_ref_log_delta;
        TF temp_ref_delta;

        Array<TF,1> press_ref, press_ref_log, temp_ref;

        Array<std::string,1> gas_names;

        Array<TF,3> vmr_ref;

        Array<int,2> flavor;
        Array<int,2> gpoint_flavor;

        Array<TF,4> kmajor;

        Array<TF,3> kminor_lower;
        Array<TF,3> kminor_upper;

        Array<int,2> minor_limits_gpt_lower;
        Array<int,2> minor_limits_gpt_upper;

        Array<int,1> minor_scales_with_density_lower;
        Array<int,1> minor_scales_with_density_upper;

        Array<int,1> scale_by_complement_lower;
        Array<int,1> scale_by_complement_upper;

        Array<int,1> kminor_start_lower;
        Array<int,1> kminor_start_upper;

        Array<int,1> idx_minor_lower;
        Array<int,1> idx_minor_upper;

        Array<int,1> idx_minor_scaling_lower;
        Array<int,1> idx_minor_scaling_upper;

        Array<int,1> is_key;

        Array<TF,1> solar_src;
        Array<TF,4> krayl;

        int get_ngas() const { return this->gas_names.dim(1); }

        void init_abs_coeffs(
                const Gas_concs<TF>& available_gases,
                const Array<std::string,1>& gas_names,
                const Array<int,3>& key_species,
                const Array<int,2>& band2gpt,
                const Array<TF,2>& band_lims_wavenum,
                const Array<TF,1>& press_ref,
                const Array<TF,1>& temp_ref,
                const TF press_ref_trop,
                const TF temp_ref_p,
                const TF temp_ref_t,
                const Array<TF,3>& vmr_ref,
                const Array<TF,4>& kmajor,
                const Array<TF,3>& kminor_lower,
                const Array<TF,3>& kminor_upper,
                const Array<std::string,1>& gas_minor,
                const Array<std::string,1>& identifier_minor,
                const Array<std::string,1>& minor_gases_lower,
                const Array<std::string,1>& minor_gases_upper,
                const Array<int,2>& minor_limits_gpt_lower,
                const Array<int,2>& minor_limits_gpt_upper,
                const Array<int,1>& minor_scales_with_density_lower,
                const Array<int,1>& minor_scales_with_density_upper,
                const Array<std::string,1>& scaling_gas_lower,
                const Array<std::string,1>& scaling_gas_upper,
                const Array<int,1>& scale_by_complement_lower,
                const Array<int,1>& scale_by_complement_upper,
                const Array<int,1>& kminor_start_lower,
                const Array<int,1>& kminor_start_upper,
                const Array<TF,3>& rayl_lower,
                const Array<TF,3>& rayl_upper);

        void compute_gas_taus(
                const int ncol, const int nlay, const int ngpt, const int nband,
                const Array<TF,2>& play,
                const Array<TF,2>& plev,
                const Array<TF,2>& tlay,
                const Gas_concs<TF>& gas_desc,
                std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                Array<int,2>& jtemp, Array<int,2>& jpress,
                Array<int,4>& jeta,
                Array<int,2>& tropo,
                Array<TF,6>& fmajor,
                const Array<TF,2>& col_dry);
};
#endif
