#ifndef RTE_LW_H
#define RTE_LW_H

#include <memory>

// Forward declarations.
template<typename, int> class Array;
template<typename> class Optical_props_arry;
template<typename> class Source_func_lw;
template<typename> class Fluxes_broadband;

template<typename TF>
class Rte_lw
{
    public:
        static void rte_lw(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const int top_at_1,
                const Source_func_lw<TF>& sources,
                const Array<TF,2>& sfc_emis,
                // CvH: RRTMGP has this: std::unique_ptr<Fluxes<TF>>& fluxes,
                std::unique_ptr<Fluxes_broadband<TF>>& fluxes,
                const int n_gauss_angles);

        static void expand_and_transpose(
                const std::unique_ptr<Optical_props_arry<TF>>& ops,
                const Array<TF,2> arr_in,
                Array<TF,2>& arr_out);
};
#endif
