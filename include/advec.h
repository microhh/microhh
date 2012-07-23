#ifndef ADVEC
#define ADVEC

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "advec_g2.h"
#include "advec_g2i4.h"
#include "advec_g42.h"

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

    cadvec_g2   *advec_g2;
    cadvec_g2i4 *advec_g2i4;
    cadvec_g42  *advec_g42;

    int iadvec;
};
#endif
