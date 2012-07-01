#ifndef BOUNDARY
#define BOUNDARY

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cboundary
{
  public:
    cboundary(cgrid *, cfields *, cmpi *);
    ~cboundary();

    int readinifile(cinput *);
    int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int setgcbot(double *, int);
    int setgctop(double *, int);

    int bcbotmom;
    int bctopmom;

    int bcbotscal;
    int bctopscal;
};
#endif
