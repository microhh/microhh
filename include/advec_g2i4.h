#ifndef ADVEC_G2I4
#define ADVEC_G2I4

#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "mpiinterface.h"

/**
 * Derived class for 2nd order advection scheme with 4th order interpolation.
 */
class cadvec_g2i4 : public cadvec
{
  public:
    cadvec_g2i4(cgrid *, cfields *, cmpi *); ///< Constructor of the advection class.
    ~cadvec_g2i4();                          ///< Destructor of the advection class.

    unsigned long gettimelim(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.
    double getcfl(double);                               ///< Get the CFL number.
    int exec();                                          ///< Execute the advection scheme.

  private:
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.
    cmpi    *mpi;    ///< Pointer to mpi class.

    double calccfl(double *, double *, double *, double *, double);         ///< Calculate the CFL number.
    int advecu(double *, double *, double *, double *, double *);           ///< Calculate longitudinal velocity advection.
    int advecv(double *, double *, double *, double *, double *);           ///< Calculate latitudinal velocity advection.
    int advecw(double *, double *, double *, double *, double *);           ///< Calculate vertical velocity advection.
    int advecs(double *, double *, double *, double *, double *, double *); ///< Calculate scalar advection.

    inline double interp2(const double, const double);                                ///< 2nd order interpolation.
    inline double interp4(const double, const double, const double, const double);    ///< 4th order interpolation.
    inline double interp4bot(const double, const double, const double, const double); ///< 4th order interpolation (bottom boundary).
    inline double interp4top(const double, const double, const double, const double); ///< 4th order interpolation (top boundary).
};
#endif
