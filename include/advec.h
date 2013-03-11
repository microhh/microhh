#ifndef ADVEC
#define ADVEC

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cadvec
{
  public:
    cadvec(cgrid *, cfields *, cmpi *);
    ~cadvec();

    int readinifile(cinput *);

    double getcfl(double);
    int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
