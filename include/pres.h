#ifndef PRES
#define PRES

#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cpres
{
  public:
    cpres(cgrid *, cfields *, cmpi *);
    ~cpres();

    int readinifile(cinput *);

    int init();
    int save();
    int load();
    int exec(double);
    double check();

  private:
    // variables
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    bool allocated;
    int ipres;

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
    int pres_2nd_solve(double *, double *, double *,
                       double *, double *, 
                       double *, double *);
    int pres_2nd_out(double *, double *, double *,
                     double *, double *);
    int tdma(double *, double *, double *, double *, 
             double *, double *);
    double calcdivergence(double *, double *, double *, double *);
};
#endif
