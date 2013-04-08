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
    virtual unsigned long gettimelim(unsigned long, double);
    virtual int exec();

    double cflmax;
  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

};
#endif
