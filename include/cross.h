#ifndef CROSS
#define CROSS

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class ccross
{
  public:
    ccross(cgrid *, cfields *, cmpi *);
    ~ccross();

    int readinifile(cinput *);
    int init();
    int exec(int);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int  icross;
    int  jxz;

    int crosssimple(double *, double *, std::string, int);
    int crosslngrad(double *, double *, double *, double *, std::string, int);
};
#endif

