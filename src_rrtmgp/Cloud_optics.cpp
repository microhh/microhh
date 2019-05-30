#include "Cloud_optics.h"

template<typename TF>
Cloud_optics<TF>::Cloud_optics(
        const Array<TF,2>& band_lims_wvn,
        const double radliq_lwr, const double radliq_upr, const double radliq_fac,
        const double radice_lwr, const double radice_upr, const double radice_fac,
        const Array<TF,2>& lut_extliq, const Array<TF,2>& lut_ssaliq, const Array<TF,2>& lut_asyliq,
        const Array<TF,3>& lut_extice, const Array<TF,3>& lut_ssaice, const Array<TF,3>& lut_asyice)
{}
