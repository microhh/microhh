#include "Cloud_optics.h"

template<typename TF>
Cloud_optics<TF>::Cloud_optics(
        const Array<TF,2>& band_lims_wvn,
        const TF radliq_lwr, const TF radliq_upr, const TF radliq_fac,
        const TF radice_lwr, const TF radice_upr, const TF radice_fac,
        const Array<TF,2>& lut_extliq, const Array<TF,2>& lut_ssaliq, const Array<TF,2>& lut_asyliq,
        const Array<TF,3>& lut_extice, const Array<TF,3>& lut_ssaice, const Array<TF,3>& lut_asyice) :
    Optical_props<TF>(band_lims_wvn)
{
    const int nsize_liq = lut_extliq.dim(1);
    const int nsize_ice = lut_extice.dim(1);

    this->liq_nsteps = nsize_liq;
    this->ice_nsteps = nsize_ice;
    this->liq_step_size = (radliq_upr - radliq_lwr) / (nsize_liq - TF(1.));
    this->ice_step_size = (radice_upr - radice_lwr) / (nsize_ice - TF(1.));

    // Load LUT constants.
    this->radliq_lwr = radliq_lwr;
    this->radliq_upr = radliq_upr;
    this->radice_lwr = radice_lwr;
    this->radice_upr = radice_upr;

    // Load LUT coefficients.
    this->lut_extliq = lut_extliq;
    this->lut_ssaliq = lut_ssaliq;
    this->lut_asyliq = lut_asyliq;
    this->lut_extice = lut_extice;
    this->lut_ssaice = lut_ssaice;
    this->lut_asyice = lut_asyice;
}

template<typename TF>
void compute_from_table(
        const int ncol, const int nlay, const int nbnd,
        const Array<int,2>& mask, const Array<TF,2>& size, const int nsteps,
        const TF step_size, const TF offset, const Array<TF,2>& table,
        Array<TF,3>& out)
{
    for (int ilay=1; ilay<=nlay; ++ilay)
        for (int icol=1; icol<=ncol; ++icol)
        {
            if (mask({icol, ilay}))
            {
                const int index = std::min(
                        static_cast<int>((size({icol, ilay}) - offset) / step_size) + 1, nsteps-1);
                const TF fint = (size({icol, ilay}) - offset) / step_size - (index-1);

                for (int ibnd=1; ibnd<=nbnd; ++ibnd)
                    out({icol, ilay, ibnd}) =
                            table({index, ibnd}) + fint * (table({index+1, ibnd}) - table({index, ibnd}));
            }
            else
            {
                for (int ibnd=1; ibnd<=nbnd; ++ibnd)
                    out({icol, ilay, ibnd}) = TF(0.);
            }
        }
}

// Two-stream variant of cloud optics.
template<typename TF>
void Cloud_optics<TF>::cloud_optics(
        const Array<int,2>& liqmsk, const Array<int,2>& icemsk,
        const Array<TF,2>& clwp, const Array<TF,2>& ciwp,
        const Array<TF,2>& reliq, const Array<TF,2>& reice,
        Optical_props_2str<TF>& optical_props)
{
    const int ncol = clwp.dim(1);
    const int nlay = clwp.dim(2);
    const int nbnd = this->get_nband();

    Optical_props_2str<TF> clouds_liq(ncol, nlay, optical_props);

    // Liquid water.
    compute_from_table(
            ncol, nlay, nbnd, liqmsk, reliq,
            this->liq_nsteps, this->liq_step_size,
            this->radliq_lwr, this->lut_extliq,
            clouds_liq.get_tau());

    compute_from_table(
            ncol, nlay, nbnd, liqmsk, reliq,
            this->liq_nsteps, this->liq_step_size,
            this->radliq_lwr, this->lut_ssaliq,
            clouds_liq.get_ssa());

    compute_from_table(
            ncol, nlay, nbnd, liqmsk, reliq,
            this->liq_nsteps, this->liq_step_size,
            this->radliq_lwr, this->lut_asyliq,
            clouds_liq.get_g());

    for (int ibnd=1; ibnd<=nbnd; ++ibnd)
        for (int ilay=1; ilay<=nlay; ++ilay)
            #pragma ivdep
            for (int icol=1; icol<=ncol; ++icol)
                clouds_liq.get_tau()({icol, ilay, ibnd}) *= clwp({icol, ilay});

    // Add ice as soon as ice microphysics are added.

    // Process the calculated optical properties.
    for (int ibnd=1; ibnd<=nbnd; ++ibnd)
        for (int ilay=1; ilay<=nlay; ++ilay)
            #pragma ivdep
            for (int icol=1; icol<=ncol; ++icol)
            {
                optical_props.get_tau()({icol, ilay, ibnd}) = clouds_liq.get_tau()({icol, ilay, ibnd});
                optical_props.get_ssa()({icol, ilay, ibnd}) = clouds_liq.get_ssa()({icol, ilay, ibnd});
                optical_props.get_g  ()({icol, ilay, ibnd}) = clouds_liq.get_g  ()({icol, ilay, ibnd});
            }
}

// 1scl variant of cloud optics.
template<typename TF>
void Cloud_optics<TF>::cloud_optics(
        const Array<int,2>& liqmsk, const Array<int,2>& icemsk,
        const Array<TF,2>& clwp, const Array<TF,2>& ciwp,
        const Array<TF,2>& reliq, const Array<TF,2>& reice,
        Optical_props_1scl<TF>& optical_props)
{
    const int ncol = clwp.dim(1);
    const int nlay = clwp.dim(2);
    const int nbnd = this->get_nband();

    Optical_props_1scl<TF> clouds_liq(ncol, nlay, optical_props);

    // Liquid water.
    compute_from_table(
            ncol, nlay, nbnd, liqmsk, reliq,
            this->liq_nsteps, this->liq_step_size,
            this->radliq_lwr, this->lut_extliq,
            clouds_liq.get_tau());

    for (int ibnd=1; ibnd<=nbnd; ++ibnd)
        for (int ilay=1; ilay<=nlay; ++ilay)
            #pragma ivdep
            for (int icol=1; icol<=ncol; ++icol)
                clouds_liq.get_tau()({icol, ilay, ibnd}) *= clwp({icol, ilay});

    // Add ice as soon as ice microphysics are added.

    // Process the calculated optical properties.
    for (int ibnd=1; ibnd<=nbnd; ++ibnd)
        for (int ilay=1; ilay<=nlay; ++ilay)
            #pragma ivdep
            for (int icol=1; icol<=ncol; ++icol)
                optical_props.get_tau()({icol, ilay, ibnd}) = clouds_liq.get_tau()({icol, ilay, ibnd});
}

#ifdef FLOAT_SINGLE
template class Cloud_optics<float>;
#else
template class Cloud_optics<double>;
#endif
