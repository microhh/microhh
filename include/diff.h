#ifndef DIFF
#define DIFF

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cdiff
{
  public:
    cdiff(cgrid *, cfields *, cmpi *);
    virtual ~cdiff();

    virtual int readinifile(cinput *);
    virtual int setvalues();
    virtual int exec();

    virtual double getdn(double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
