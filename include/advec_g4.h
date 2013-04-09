#ifndef ADVEC_G4
#define ADVEC_G4

#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "mpiinterface.h"

/**
 * Derived class for 4th order advection scheme.
 */
class cadvec_g4 : public cadvec
{
  public:
    cadvec_g4(cgrid *, cfields *, cmpi *); ///< Constructor of the advection class.
    ~cadvec_g4();                          ///< Destructor of the advection class.

    unsigned long gettimelim(long unsigned int idt, double ifactor);
    double getcfl(double);
    int exec();

  private:
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.
    cmpi    *mpi;    ///< Pointer to mpi class.

    double calccfl(double *, double *, double *, double *, double);         ///< Calculate the CFL number.
    int advecu(double *, double *, double *, double *, double *);           ///< Calculate longitudinal velocity advection.
    int advecv(double *, double *, double *, double *, double *);           ///< Calculate latitudinal velocity advection.
    int advecw(double *, double *, double *, double *, double *);           ///< Calculate vertical velocity advection.
    int advecs(double *, double *, double *, double *, double *, double *); ///< Calculate scalar advection.
};
#endif
