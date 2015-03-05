#ifndef THERMO_MOIST_FUNCTIONS
#define THERMO_MOIST_FUNCTIONS

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

#include "constants.h"

namespace Thermo_moist_functions
{
    using namespace constants;

    // INLINE FUNCTIONS
    CUDA_MACRO inline double buoyancy(const double exn, const double thl, const double qt, const double ql, const double thvref)
    {
        return grav * ((thl + Lv*ql/(cp*exn)) * (1. - (1. - Rv/Rd)*qt - Rv/Rd*ql) - thvref) / thvref;
    }

    CUDA_MACRO inline double buoyancy_no_ql(const double thl, const double qt, const double thvref)
    {
        return grav * (thl * (1. - (1. - Rv/Rd)*qt) - thvref) / thvref;
    }

    CUDA_MACRO inline double buoyancy_flux_no_ql(const double thl, const double thlflux, const double qt, const double qtflux, const double thvref)
    {
        return grav/thvref * (thlflux * (1. - (1.-Rv/Rd)*qt) - (1.-Rv/Rd)*thl*qtflux);
    }

    CUDA_MACRO inline double esat(const double T)
    {
        const double x=std::max(-80.,T-T0);
        return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
    }

    CUDA_MACRO inline double qsat(const double p, const double T)
    {
        return ep*esat(T)/(p-(1-ep)*esat(T));
    }

    CUDA_MACRO inline double sat_adjust(const double thl, const double qt, const double p, const double exn)
    {
        int niter = 0;
        int nitermax = 30;
        double ql, tl, tnr_old = 1.e9, tnr, qs=0;
        tl = thl * exn;
        tnr = tl;
        while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5 && niter < nitermax)
        {
            ++niter;
            tnr_old = tnr;
            qs = qsat(p,tnr);
            tnr = tnr - (tnr+(Lv/cp)*qs-tl-(Lv/cp)*qt)/(1+(std::pow(Lv,2)*qs)/ (Rv*cp*std::pow(tnr,2)));
        }

        if (niter == nitermax)
        {  
            printf("Saturation adjustment not converged!! [thl=%f K, qt=%f kg/kg, p=%f p]\n",thl,qt,p);
            throw 1;
        }  

        ql = std::max(0.,qt - qs);
        return ql;
    }

    CUDA_MACRO inline double exner(const double p)
    {
        return pow((p/p0),(Rd/cp));
    }
}
#endif
