#ifndef DIFF_G42
#define DIFF_G42

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cdiff_g42
{
  public:
    cdiff_g42(cgrid *, cfields *, cmpi *);
    ~cdiff_g42();

    int diffc(double *, double *, double *, double *, double);
    int diffw(double *, double *, double *, double *, double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    inline double divgrad4(const double, const double, const double, const double,
                           const double, const double, const double, const double);
};
#endif
