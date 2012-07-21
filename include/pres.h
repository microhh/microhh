#ifndef PRES
#define PRES

#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "pres_g2.h"

class cpres
{
  public:
    cpres(cgrid *, cfields *, cmpi *);
    ~cpres();

    int readinifile(cinput *);

    int init();
    int load();
    int save();

    int exec(double);

    double check();

  private:
    // variables
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int ipres;

    cpres_g2 *pres_g2;
};
#endif
