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

#ifndef SOURCE_FUNCTIONS_H
#define SOURCE_FUNCTIONS_H

template<typename, int> class Array;
template<typename> class Optical_props;
template<typename> class Source_func_lw;

template<typename TF>
class Source_func_lw : public Optical_props<TF>
{
    public:
        Source_func_lw(
                const int n_col,
                const int n_lay,
                const Optical_props<TF>& optical_props);

        void set_subset(
                const Source_func_lw<TF>& sources_sub,
                const int col_s, const int col_e);

        void get_subset(
                const Source_func_lw<TF>& sources_sub,
                const int col_s, const int col_e);

        Array<TF,2>& get_sfc_source()     { return sfc_source;     }
        Array<TF,2>& get_sfc_source_jac() { return sfc_source_jac; }
        Array<TF,3>& get_lay_source()     { return lay_source;     }
        Array<TF,3>& get_lev_source_inc() { return lev_source_inc; }
        Array<TF,3>& get_lev_source_dec() { return lev_source_dec; }

        const Array<TF,2>& get_sfc_source()     const { return sfc_source;     }
        const Array<TF,2>& get_sfc_source_jac() const { return sfc_source_jac; }
        const Array<TF,3>& get_lay_source()     const { return lay_source;     }
        const Array<TF,3>& get_lev_source_inc() const { return lev_source_inc; }
        const Array<TF,3>& get_lev_source_dec() const { return lev_source_dec; }

    private:
        Array<TF,2> sfc_source;
        Array<TF,2> sfc_source_jac;
        Array<TF,3> lay_source;
        Array<TF,3> lev_source_inc;
        Array<TF,3> lev_source_dec;
};
#endif
