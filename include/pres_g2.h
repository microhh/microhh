#ifndef PRES_G2
#define PRES_G2

#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "mpiinterface.h"

class cpres_g2 : public cpres
{
  public:
    cpres_g2(cgrid *, cfields *, cmpi *);
    ~cpres_g2();

    int init();
    int setvalues();
    int exec(double);
    double check();

  private:
    // variables
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    bool allocated;

    double *bmati, *bmatj;
    double *a, *c;
    double *work2d;

    int pres_in(double *, 
                double *, double *, double *,
                double *, double *, double *,
                double *, double);
    int pres_solve(double *, double *, double *, double *,
                   double *, double *, 
                   double *, double *);
    int pres_out(double *, double *, double *,
                 double *, double *);
    double calcdivergence(double *, double *, double *, double *);

    // functions
    int tdma(double *, double *, double *, double *, 
             double *, double *);
};
#endif
