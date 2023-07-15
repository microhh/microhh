#ifndef ADVEC_MONOTONIC_H
#define ADVEC_MONOTONIC_H

#include <cmath>

namespace Advec_monotonic
{
    // Implementation flux limiter according to Koren, 1993.
    template<typename TF>
    inline TF flux_lim(const TF u, const TF sm2, const TF sm1, const TF sp1, const TF sp2)
    {
        const TF eps = std::numeric_limits<TF>::epsilon();

        if (u >= TF(0.))
        {
            const TF denom = copysign(1, sm1-sm2) * std::max(std::abs(sm1-sm2), eps);
            const TF two_r = TF(2.) * (sp1-sm1) / denom;
            const TF phi = std::max(
                    TF(0.),
                    std::min( two_r, std::min( TF(1./3.)*(TF(1.)+two_r), TF(2.)) ) );
            return u*(sm1 + TF(0.5)*phi*(sm1 - sm2));
        }
        else
        {
            const TF denom = copysign(1, sp1-sp2) * std::max(std::abs(sp1-sp2), eps);
            const TF two_r = TF(2.) * (sm1-sp1) / denom;
            const TF phi = std::max(
                    TF(0.),
                    std::min( two_r, std::min( TF(1./3.)*(TF(1.)+two_r), TF(2.)) ) );
            return u*(sp1 + TF(0.5)*phi*(sp1 - sp2));
        }
    }

    // Implementation flux limiter according to Koren, 1993.
    template<typename TF>
    inline TF flux_lim_bot(const TF u, const TF sm2, const TF sm1, const TF sp1, const TF sp2)
    {
        const TF eps = std::numeric_limits<TF>::epsilon();

        if (u >= TF(0.))
        {
            return u*sm1;
        }
        else
        {
            const TF denom = copysign(1, sp1-sp2) * std::max(std::abs(sp1-sp2), eps);
            const TF two_r = TF(2.) * (sm1-sp1) / denom;
            const TF phi = std::max(
                    TF(0.),
                    std::min( two_r, std::min( TF(1./3.)*(TF(1.)+two_r), TF(2.)) ) );
            return u*(sp1 + TF(0.5)*phi*(sp1 - sp2));
        }
    }

    // Implementation flux limiter according to Koren, 1993.
    template<typename TF>
    inline TF flux_lim_top(const TF u, const TF sm2, const TF sm1, const TF sp1, const TF sp2)
    {
        const TF eps = std::numeric_limits<TF>::epsilon();

        if (u >= TF(0.))
        {
            const TF denom = copysign(1, sm1-sm2) * std::max(std::abs(sm1-sm2), eps);
            const TF two_r = TF(2.) * (sp1-sm1) / denom;
            const TF phi = std::max(
                    TF(0.),
                    std::min( two_r, std::min( TF(1./3.)*(TF(1.)+two_r), TF(2.)) ) );
            return u*(sm1 + TF(0.5)*phi*(sm1 - sm2));
        }
        else
        {
            return u*sp1;
        }
    }

    template<typename TF>
    void advec_s_lim(
            TF* const restrict st, const TF* const restrict s,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        for (int k=kstart+2; k<kend-2; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    st[ijk] +=
                             - ( flux_lim(u[ijk+ii1], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                               - flux_lim(u[ijk    ], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                             - ( flux_lim(v[ijk+jj1], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                               - flux_lim(v[ijk    ], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                             - ( rhorefh[k+1] * flux_lim(w[ijk+kk1], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                               - rhorefh[k  ] * flux_lim(w[ijk    ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
                }

        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( flux_lim(u[ijk+ii1], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - flux_lim(u[ijk    ], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( flux_lim(v[ijk+jj1], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - flux_lim(v[ijk    ], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         // No flux through bottom wall.
                         - ( rhorefh[k+1] * flux_lim_bot(w[ijk+kk1], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])) / rhoref[k] * dzi[k];
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( flux_lim(u[ijk+ii1], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - flux_lim(u[ijk    ], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( flux_lim(v[ijk+jj1], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - flux_lim(v[ijk    ], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         - ( rhorefh[k+1] * flux_lim    (w[ijk+kk1], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                           - rhorefh[k  ] * flux_lim_bot(w[ijk    ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        k = kend-2;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( flux_lim(u[ijk+ii1], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - flux_lim(u[ijk    ], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( flux_lim(v[ijk+jj1], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - flux_lim(v[ijk    ], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         // No flux through bottom wall.
                         - ( rhorefh[k+1] * flux_lim_top(w[ijk+kk1], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                           - rhorefh[k  ] * flux_lim    (w[ijk    ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( flux_lim(u[ijk+ii1], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - flux_lim(u[ijk    ], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( flux_lim(v[ijk+jj1], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - flux_lim(v[ijk    ], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         - ( // No flux through boundary
                           - rhorefh[k  ] * flux_lim_top(w[ijk    ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
    }

    template<typename TF>
    void advec_flux_s_lim(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] = TF(0.);
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] = flux_lim_bot(w[ijk], s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = flux_lim(w[ijk], s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] = flux_lim_top(w[ijk], s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
            }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] = TF(0.);
            }
    }
}
#endif
