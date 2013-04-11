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

    virtual unsigned long gettimelim(unsigned long, double);
    virtual double getdn(double);

    double dnmax;
  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
