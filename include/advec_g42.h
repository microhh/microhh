#ifndef ADVEC_G42
#define ADVEC_G42

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cadvec_g42
{
  public:
    cadvec_g42(cgrid *, cfields *, cmpi *);
    ~cadvec_g42();

    double calccfl(double *, double *, double *, double *, double);
    int advecu(double *, double *, double *, double *, double *);
    int advecv(double *, double *, double *, double *, double *);
    int advecw(double *, double *, double *, double *, double *);
    int advecs(double *, double *, double *, double *, double *, double *);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    inline double interp2(const double, const double);
};
#endif
