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

#ifdef FLOAT_SINGLE
template class Gas_concs<float>;
#else
template class Gas_concs<double>;
#endif
