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

#ifndef FLUXES_H
#define FLUXES_H

#include <memory>
#include <stdexcept>

#include "define_bool.h"

// Forward declarations.
template<typename TF, int> class Array;
template<typename TF> class Optical_props_arry;

template<typename TF>
class Fluxes
{
    public:
        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const BOOL_TYPE top_at_1) = 0;

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const Array<TF,3>& gpt_flux_dn_dir,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const BOOL_TYPE top_at_1) = 0;
};

template<typename TF>
class Fluxes_broadband : public Fluxes<TF>
{
    public:
        Fluxes_broadband(const int ncol, const int nlev);
        virtual ~Fluxes_broadband() {};

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const BOOL_TYPE top_at_1);

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const Array<TF,3>& gpt_flux_dn_dir,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const BOOL_TYPE top_at_1);

        Array<TF,2>& get_flux_up    () { return flux_up;     }
        Array<TF,2>& get_flux_dn    () { return flux_dn;     }
        Array<TF,2>& get_flux_dn_dir() { return flux_dn_dir; }
        Array<TF,2>& get_flux_net   () { return flux_net;    }

        virtual Array<TF,3>& get_bnd_flux_up    () { throw std::runtime_error("Band fluxes are not available"); }
        virtual Array<TF,3>& get_bnd_flux_dn    () { throw std::runtime_error("Band fluxes are not available"); }
        virtual Array<TF,3>& get_bnd_flux_dn_dir() { throw std::runtime_error("Band fluxes are not available"); }
        virtual Array<TF,3>& get_bnd_flux_net   () { throw std::runtime_error("Band fluxes are not available"); }

    private:
        Array<TF,2> flux_up;
        Array<TF,2> flux_dn;
        Array<TF,2> flux_dn_dir;
        Array<TF,2> flux_net;
};

template<typename TF>
class Fluxes_byband : public Fluxes_broadband<TF>
{
    public:
        Fluxes_byband(const int ncol, const int nlev, const int nbnd);
        virtual ~Fluxes_byband() {};

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const BOOL_TYPE top_at_1);

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const Array<TF,3>& gpt_flux_dn_dir,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const BOOL_TYPE top_at_1);

        Array<TF,3>& get_bnd_flux_up    () { return bnd_flux_up;     }
        Array<TF,3>& get_bnd_flux_dn    () { return bnd_flux_dn;     }
        Array<TF,3>& get_bnd_flux_dn_dir() { return bnd_flux_dn_dir; }
        Array<TF,3>& get_bnd_flux_net   () { return bnd_flux_net;    }

    private:
        Array<TF,3> bnd_flux_up;
        Array<TF,3> bnd_flux_dn;
        Array<TF,3> bnd_flux_dn_dir;
        Array<TF,3> bnd_flux_net;
};
#endif
