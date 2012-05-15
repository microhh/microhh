#ifndef DIFF
#define DIFF

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cdiff
{
  public:
    cdiff(cgrid *, cfields *, cmpi *);
    ~cdiff();

    int init();
    int exec();

    double getdn(double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int diffc_2nd(double *, double *, double *, double *, double);
    int diffw_2nd(double *, double *, double *, double *, double);

    double dnmul;
};
#endif
