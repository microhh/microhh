#ifndef CLOUD_OPTICS_H
#define CLOUD_OPTICS_H

#include "Array.h"
#include "Optical_props.h"

// Forward declarations.
template<typename TF> class Optical_props;

template<typename TF>
class Cloud_optics : public Optical_props<TF>
{
    public:
        Cloud_optics(
                const Array<TF,2>& band_lims_wvn,
                const TF radliq_lwr, const TF radliq_upr, const TF radliq_fac,
                const TF radice_lwr, const TF radice_upr, const TF radice_fac,
                const Array<TF,2>& lut_extliq, const Array<TF,2>& lut_ssaliq, const Array<TF,2>& lut_asyliq,
                const Array<TF,3>& lut_extice, const Array<TF,3>& lut_ssaice, const Array<TF,3>& lut_asyice);
        void cloud_optics() {}

    private:
        int liq_nsteps;
        int ice_nsteps;
        TF liq_step_size;
        TF ice_step_size;

        // Lookup table constants.
        TF radliq_lwr;
        TF radliq_upr;
        TF radice_lwr;
        TF radice_upr;

        // Lookup table coefficients.
        Array<TF,2> lut_extliq;
        Array<TF,2> lut_ssaliq;
        Array<TF,2> lut_asyliq;
        Array<TF,3> lut_extice;
        Array<TF,3> lut_ssaice;
        Array<TF,3> lut_asyice;
};
#endif
