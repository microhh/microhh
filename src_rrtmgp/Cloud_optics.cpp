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
void Cloud_optics<TF>::cloud_optics(
        const int ncol, const int nlay, const int nbnd, const int nrghice,
        const Array<int,2>& liqmsk, const Array<int,2>& icemsk,
        const Array<TF,2>& clwp, const Array<int,2>& ciwp,
        const Array<TF,2>& reliq, const Array<int,2>& reice,
        std::unique_ptr<Optical_props_arry<TF>>& optical_props)
{
}

#ifdef FLOAT_SINGLE
template class Cloud_optics<float>;
#else
template class Cloud_optics<double>;
#endif
