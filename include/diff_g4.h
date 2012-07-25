#ifndef DIFF_G4
#define DIFF_G4

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cdiff_g4
{
  public:
    cdiff_g4(cgrid *, cfields *, cmpi *);
    ~cdiff_g4();

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
