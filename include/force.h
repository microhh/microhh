#ifndef FORCE
#define FORCE

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cforce
{
  public:
    cforce(cgrid *, cfields *, cmpi *);
    ~cforce();
    int readinifile(cinput *);
    int exec(double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int iforce;
    double uflow;

    int flux(double *, double *, double *, double);
};
#endif
