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

#ifndef RTE_SW_H
#define RTE_SW_H

#include <memory>

// Forward declarations.
template<typename, int> class Array;
template<typename> class Optical_props_arry;
template<typename> class Fluxes_broadband;

template<typename TF>
class Rte_sw
{
    public:
        static void rte_sw(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const BOOL_TYPE top_at_1,
                const Array<TF,1>& mu0,
                const Array<TF,2>& inc_flux_dir,
                const Array<TF,2>& sfc_alb_dir,
                const Array<TF,2>& sfc_alb_dif,
                const Array<TF,2>& inc_flux_dif,
                Array<TF,3>& gpt_flux_up,
                Array<TF,3>& gpt_flux_dn,
                Array<TF,3>& gpt_flux_dir);

        static void expand_and_transpose(
                const std::unique_ptr<Optical_props_arry<TF>>& ops,
                const Array<TF,2> arr_in,
                Array<TF,2>& arr_out);
};
#endif
