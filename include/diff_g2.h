#ifndef DIFF_G2
#define DIFF_G2

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cdiff_g2
{
  public:
    cdiff_g2(cgrid *, cfields *, cmpi *);
    ~cdiff_g2();

    int diffc(double *, double *, double *, double *, double);
    int diffw(double *, double *, double *, double *, double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
