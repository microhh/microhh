#ifndef PRES
#define PRES

#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cpres
{
  public:
    cpres(cgrid *, cfields *, cmpi *);
    virtual ~cpres();

    virtual int readinifile(cinput *);

    virtual int init();
    virtual int setvalues();
    virtual int exec(double);
    virtual double check();

  private:
    // variables
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
