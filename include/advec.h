#ifndef ADVEC
#define ADVEC

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cadvec
{
  public:
    cadvec(cgrid *, cfields *, cmpi *);
    virtual ~cadvec();

    virtual int readinifile(cinput *);

    virtual double getcfl(double);
    virtual int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
