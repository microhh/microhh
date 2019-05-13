#include "Source_functions.h"
#include "Array.h"
#include "Optical_props.h"

template<typename TF>
Source_func_lw<TF>::Source_func_lw(
        const int n_col,
        const int n_lay,
        const Optical_props<TF>& optical_props) :
    Optical_props<TF>(optical_props),
    sfc_source({n_col, optical_props.get_ngpt()}),
    lay_source({n_col, n_lay, optical_props.get_ngpt()}),
    lev_source_inc({n_col, n_lay, optical_props.get_ngpt()}),
    lev_source_dec({n_col, n_lay, optical_props.get_ngpt()})
{}

template<typename TF>
void Source_func_lw<TF>::set_subset(
        const Source_func_lw<TF>& sources_sub,
        const int col_s, const int col_e)
{
    for (int igpt=1; igpt<=lay_source.dim(3); ++igpt)
        for (int icol=col_s; icol<=col_e; ++icol)
            sfc_source({icol, igpt}) = sources_sub.get_sfc_source()({icol-col_s+1, igpt});

    for (int igpt=1; igpt<=lay_source.dim(3); ++igpt)
        for (int ilay=1; ilay<=lay_source.dim(2); ++ilay)
            for (int icol=col_s; icol<=col_e; ++icol)
            {
                lay_source    ({icol, ilay, igpt}) = sources_sub.get_lay_source()    ({icol-col_s+1, ilay, igpt});
                lev_source_inc({icol, ilay, igpt}) = sources_sub.get_lev_source_inc()({icol-col_s+1, ilay, igpt});
                lev_source_dec({icol, ilay, igpt}) = sources_sub.get_lev_source_dec()({icol-col_s+1, ilay, igpt});
            }
}

template<typename TF>
void Source_func_lw<TF>::get_subset(
        const Source_func_lw<TF>& sources_sub,
        const int col_s, const int col_e)
{
    for (int igpt=1; igpt<=lay_source.dim(3); ++igpt)
        for (int icol=col_s; icol<=col_e; ++icol)
            sfc_source({icol-col_s+1, igpt}) = sources_sub.get_sfc_source()({icol, igpt});

    for (int igpt=1; igpt<=lay_source.dim(3); ++igpt)
        for (int ilay=1; ilay<=lay_source.dim(2); ++ilay)
            for (int icol=col_s; icol<=col_e; ++icol)
            {
                lay_source    ({icol-col_s+1, ilay, igpt}) = sources_sub.get_lay_source()    ({icol, ilay, igpt});
                lev_source_inc({icol-col_s+1, ilay, igpt}) = sources_sub.get_lev_source_inc()({icol, ilay, igpt});
                lev_source_dec({icol-col_s+1, ilay, igpt}) = sources_sub.get_lev_source_dec()({icol, ilay, igpt});
            }
}

#ifdef FLOAT_SINGLE
template class Source_func_lw<float>;
#else
template class Source_func_lw<double>;
#endif
