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

    int setgcbot_2nd(double *, int);
    int setgctop_2nd(double *, int);
    int setgcbot_4th(double *, int);
    int setgctop_4th(double *, int);

    int iboundary;

    int bcbotmom;
    int bctopmom;

    int bcbotscal;
    int bctopscal;
};
#endif
