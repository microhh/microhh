#ifndef ADVEC_G4M
#define ADVEC_G4M

#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "mpiinterface.h"

/**
 * Derived class for fully conservative 4th order advection scheme.
 * Fully mass, momentum and energy conserving advection scheme based on the paper
 * of Morinishi et al., (1998).
 */
class cadvec_g4m : public cadvec
{
  public:
    cadvec_g4m(cgrid *, cfields *, cmpi *); ///< Constructor of the advection class.
    ~cadvec_g4m();                          ///< Destructor of the advection class.

    unsigned long gettimelim(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.
    double getcfl(double);                               ///< Get the CFL number.
    int exec();                                          ///< Execute the advection scheme.

  private:
    double calccfl(double *, double *, double *, double *, double);         ///< Calculate the CFL number.
    int advecu(double *, double *, double *, double *, double *);           ///< Calculate longitudinal velocity advection.
    int advecv(double *, double *, double *, double *, double *);           ///< Calculate latitudinal velocity advection.
    int advecw(double *, double *, double *, double *, double *);           ///< Calculate vertical velocity advection.
    int advecs(double *, double *, double *, double *, double *, double *); ///< Calculate scalar advection.

    inline double interp2(const double, const double);                                           ///< 2nd order interpolation.
    inline double interp4(const double, const double, const double, const double);               ///< 4th order interpolation.
    inline double grad4  (const double, const double, const double, const double, const double); ///< 4th order gradient.
    inline double grad4x (const double, const double, const double, const double);               ///< 4th order gradient (only numerator).

    inline double interp4biasbot(const double, const double, const double, const double); ///< 4th order interpolation (bottom boundary).
    inline double interp4biastop(const double, const double, const double, const double); ///< 4th order interpolation (top boundary).
    inline double grad4xbiasbot (const double, const double, const double, const double); ///< 4th order interpolation (bottom boundary, only numerator).
    inline double grad4xbiastop (const double, const double, const double, const double); ///< 4th order interpolation (top boundary, only numerator).
};
#endif
