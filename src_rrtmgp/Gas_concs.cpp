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

#include "Gas_concs.h"
#include "Array.h"

template<typename TF>
Gas_concs<TF>::Gas_concs(const Gas_concs& gas_concs_ref, const int start, const int size)
{
    const int end = start + size - 1;
    for (auto& g : gas_concs_ref.gas_concs_map)
    {
        if (g.second.dim(1) == 1)
            this->gas_concs_map.emplace(g.first, g.second);
        else
        {
            Array<TF,2> gas_conc_subset = g.second.subset({{ {start, end}, {1, g.second.dim(2)} }});
            this->gas_concs_map.emplace(g.first, gas_conc_subset);
        }
    }
}

// Insert new gas into the map or update the value.
template<typename TF>
void Gas_concs<TF>::set_vmr(const std::string& name, const TF data)
{
    Array<TF,2> data_2d({1, 1});
    data_2d({1, 1}) = data;

    if (this->exists(name))
        gas_concs_map.at(name) = data_2d;
    else
        gas_concs_map.emplace(name, std::move(data_2d));
}

// Insert new gas into the map or update the value.
template<typename TF>
void Gas_concs<TF>::set_vmr(const std::string& name, const Array<TF,1>& data)
{
    Array<TF,2> data_2d(data.v(), {1, data.dim(1)});

    if (this->exists(name))
        gas_concs_map.at(name) = data_2d;
    else
        gas_concs_map.emplace(name, std::move(data_2d));
}

// Insert new gas into the map or update the value.
template<typename TF>
void Gas_concs<TF>::set_vmr(const std::string& name, const Array<TF,2>& data_2d)
{
    if (this->exists(name))
        gas_concs_map.at(name) = data_2d;
    else
        gas_concs_map.emplace(name, data_2d);
}

// Get gas from map.
template<typename TF>
const Array<TF,2>& Gas_concs<TF>::get_vmr(const std::string& name) const
{
    return this->gas_concs_map.at(name);
}

// Check if gas exists in map.
template<typename TF>
bool Gas_concs<TF>::exists(const std::string& name) const
{ 
    return gas_concs_map.count(name) != 0;
}

#ifdef FLOAT_SINGLE_RRTMGP
template class Gas_concs<float>;
#else
template class Gas_concs<double>;
#endif
