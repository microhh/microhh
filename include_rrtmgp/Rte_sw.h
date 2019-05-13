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
                const int top_at_1,
                const Array<TF,1>& mu0,
                const Array<TF,2>& inc_flux,
                const Array<TF,2>& sfc_alb_dir,
                const Array<TF,2>& sfc_alb_dif,
                std::unique_ptr<Fluxes_broadband<TF>>& fluxes);

        static void expand_and_transpose(
                const std::unique_ptr<Optical_props_arry<TF>>& ops,
                const Array<TF,2> arr_in,
                Array<TF,2>& arr_out);
};
#endif
