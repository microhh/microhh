#include "Optical_props.h"
#include "Array.h"

template<typename TF>
Optical_props<TF>::Optical_props(
        Array<TF,2>& band_lims_wvn,
        Array<int,2>& band_lims_gpt)
{
    Array<int,2> band_lims_gpt_lcl(band_lims_gpt);
    this->band2gpt = band_lims_gpt_lcl;
    this->band_lims_wvn = band_lims_wvn;

    // Make a map between g-points and bands.
    this->gpt2band.set_dims({band_lims_gpt_lcl.max()});
    for (int iband=1; iband<=band_lims_gpt_lcl.dim(2); ++iband)
    {
        for (int i=band_lims_gpt_lcl({1,iband}); i<=band_lims_gpt_lcl({2,iband}); ++i)
            this->gpt2band({i}) = iband;
    }
}

template<typename TF>
Optical_props_1scl<TF>::Optical_props_1scl(
        const int ncol,
        const int nlay,
        const Optical_props<TF>& optical_props) :
    Optical_props_arry<TF>(optical_props),
    tau({ncol, nlay, this->get_ngpt()})
{}

template<typename TF>
void Optical_props_1scl<TF>::set_subset(
        const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
        const int col_s, const int col_e)
{
    for (int igpt=1; igpt<=tau.dim(3); ++igpt)
        for (int ilay=1; ilay<=tau.dim(2); ++ilay)
            for (int icol=col_s; icol<=col_e; ++icol)
                tau({icol, ilay, igpt}) = optical_props_sub->get_tau()({icol-col_s+1, ilay, igpt});
}

template<typename TF>
void Optical_props_1scl<TF>::get_subset(
        const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
        const int col_s, const int col_e)
{
    for (int igpt=1; igpt<=tau.dim(3); ++igpt)
        for (int ilay=1; ilay<=tau.dim(2); ++ilay)
            for (int icol=col_s; icol<=col_e; ++icol)
                tau({icol-col_s+1, ilay, igpt}) = optical_props_sub->get_tau()({icol, ilay, igpt});
}

template<typename TF>
Optical_props_2str<TF>::Optical_props_2str(
        const int ncol,
        const int nlay,
        const Optical_props<TF>& optical_props) :
    Optical_props_arry<TF>(optical_props),
    tau({ncol, nlay, this->get_ngpt()}),
    ssa({ncol, nlay, this->get_ngpt()}),
    g  ({ncol, nlay, this->get_ngpt()})
{}

template<typename TF>
void Optical_props_2str<TF>::set_subset(
        const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
        const int col_s, const int col_e)
{
    for (int igpt=1; igpt<=tau.dim(3); ++igpt)
        for (int ilay=1; ilay<=tau.dim(2); ++ilay)
            for (int icol=col_s; icol<=col_e; ++icol)
            {
                tau({icol, ilay, igpt}) = optical_props_sub->get_tau()({icol-col_s+1, ilay, igpt});
                ssa({icol, ilay, igpt}) = optical_props_sub->get_ssa()({icol-col_s+1, ilay, igpt});
                g  ({icol, ilay, igpt}) = optical_props_sub->get_g  ()({icol-col_s+1, ilay, igpt});
            }
}

template<typename TF>
void Optical_props_2str<TF>::get_subset(
        const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
        const int col_s, const int col_e)
{
    for (int igpt=1; igpt<=tau.dim(3); ++igpt)
        for (int ilay=1; ilay<=tau.dim(2); ++ilay)
            for (int icol=col_s; icol<=col_e; ++icol)
            {
                tau({icol-col_s+1, ilay, igpt}) = optical_props_sub->get_tau()({icol, ilay, igpt});
                ssa({icol-col_s+1, ilay, igpt}) = optical_props_sub->get_ssa()({icol, ilay, igpt});
                g  ({icol-col_s+1, ilay, igpt}) = optical_props_sub->get_g  ()({icol, ilay, igpt});
            }
}

#ifdef FLOAT_SINGLE
template class Optical_props<float>;
template class Optical_props_1scl<float>;
template class Optical_props_2str<float>;
#else
template class Optical_props<double>;
template class Optical_props_1scl<double>;
template class Optical_props_2str<double>;
#endif
