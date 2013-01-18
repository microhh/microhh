#ifndef DIFF_LES_G2
#define DIFF_LES_G2

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cdiff_les_g2
{
  public:
    cdiff_les_g2(cgrid *, cfields *, cmpi *);
    ~cdiff_les_g2();

    int evisc(double *, double *, double *, double *, double *, double *, double *, double);
    int diffu(double *, double *, double *, double *, double *, double *, double *);
    int diffv(double *, double *, double *, double *, double *, double *, double *);
    int diffw(double *, double *, double *, double *, double *, double *, double *);
    int diffc(double *, double *, double *, double *, double *);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
