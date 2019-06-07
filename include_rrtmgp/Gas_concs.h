#ifndef GAS_CONCS_H
#define GAS_CONCS_H

#include <map>

template<typename, int> class Array;

template<typename TF>
class Gas_concs
{
    public:
        Gas_concs() {}
        Gas_concs(const Gas_concs& gas_concs_ref, const int start, const int size);

        // Insert new gas into the map.
        void set_vmr(const std::string& name, const TF data);
        void set_vmr(const std::string& name, const Array<TF,1>& data);
        void set_vmr(const std::string& name, const Array<TF,2>& data);

        // Insert new gas into the map.
        // void get_vmr(const std::string& name, Array<TF,2>& data) const;
        const Array<TF,2>& get_vmr(const std::string& name) const;

        // Check if gas exists in map.
        bool exists(const std::string& name) const;

    private:
        std::map<std::string, Array<TF,2>> gas_concs_map;
};
#endif
