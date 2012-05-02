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

    int exec(double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
    int flux(double *, double *, double *, double);

};
#endif
