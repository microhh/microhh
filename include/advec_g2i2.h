#ifndef ADVEC_G2I2
#define ADVEC_G2I2

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cadvec
{
  public:
    cadvec(cgrid *, cfields *, cmpi *);
    ~cadvec();

    int exec();
    double getcfl(double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double calccfl(double *, double *, double *, double *, double);
    int advecu(double *, double *, double *, double *, double *);
    int advecv(double *, double *, double *, double *, double *);
    int advecw(double *, double *, double *, double *, double *);
    int advecs(double *, double *, double *, double *, double *, double *);
    inline double interp2(const double, const double);
};
#endif
