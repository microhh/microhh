/*
 * This file is part of a C++ interface to the Radiative Transfer for Energetics (RTE)
 * and Rapid Radiative Transfer Model for GCM applications Parallel (RRTMGP).
 *
 * The original code is found at https://github.com/RobertPincus/rte-rrtmgp.
 *
 * Contacts: Robert Pincus and Eli Mlawer
 * email: rrtmgp@aer.com
 *
 * Copyright 2015-2019,  Atmospheric and Environmental Research and
 * Regents of the University of Colorado.  All right reserved.
 *
 * This C++ interface can be downloaded from https://github.com/microhh/rrtmgp_cpp
 *
 * Contact: Chiel van Heerwaarden
 * email: chiel.vanheerwaarden@wur.nl
 *
 * Copyright 2019, Wageningen University & Research.
 *
 * Use and duplication is permitted under the terms of the
 * BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
 *
 */

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

        void cloud_optics(
                const Array<int,2>& liqmsk, const Array<int,2>& icemsk,
                const Array<TF,2>& clwp, const Array<TF,2>& ciwp,
                const Array<TF,2>& reliq, const Array<TF,2>& reice,
                Optical_props_1scl<TF>& optical_props);

        void cloud_optics(
                const Array<int,2>& liqmsk, const Array<int,2>& icemsk,
                const Array<TF,2>& clwp, const Array<TF,2>& ciwp,
                const Array<TF,2>& reliq, const Array<TF,2>& reice,
                Optical_props_2str<TF>& optical_props);

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
