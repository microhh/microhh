#ifndef PRES
#define PRES

#include <fftw3.h>
#include "grid.h"
#include "fields.h"

class cpres
{
  public:
    cpres(cgrid *, cfields *);
    ~cpres();

    int init();
    int exec(double);
    int divergence();

  private:
    // variables
    cgrid *grid;
    cfields *fields;

    double *fftini, *fftouti;
    double *fftinj, *fftoutj;

    fftw_plan iplanf, iplanb;
    fftw_plan jplanf, jplanb;

    double *bmati, *bmatj;

    double *a, *b, *c;
    double *work2d, *work3d;

    // functions
    int pres_2nd_init();
    int pres_2nd_in(double *, 
                    double *, double *, double *,
                    double *, double *, double *,
                    double *, double);
    int pres_2nd_solve(double *, double *);
    int pres_2nd_out(double *, double *, double *,
                     double *, double *);
    int tdma(double *, double *, double *, double *, 
             double *, double *);
    double calcdivergence(double *, double *, double *, double *);
};
#endif
