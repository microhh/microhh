#ifndef PRES_G4
#define PRES_G4

#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cpres_g4
{
  public:
    cpres_g4(cgrid *, cfields *, cmpi *);
    ~cpres_g4();

    int init();
    int load();
    int save();
    int setvalues();

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

    double *fftini, *fftouti;
    double *fftinj, *fftoutj;

  private:
    // variables
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    bool allocated;

    fftw_plan iplanf, iplanb;
    fftw_plan jplanf, jplanb;

    double *bmati, *bmatj;
    double *m0,*m1,*m2,*m3,*m4,*m5,*m6,*m7,*m8;
    double *work2d;

    // CvH remove later
    double *m0temp,*m1temp,*m2temp,*m3temp,*m4temp,*m5temp,*m6temp,*m7temp,*m8temp,*ptemp;

    // functions
    int hdma(double *, double *, double *, double *,
             double *, double *, double *, double *);

    /*
    inline double grad4(const double, const double, const double, const double, const double);

    inline double grad4x       (const double, const double, const double, const double);
    inline double grad4xbiasbot(const double, const double, const double, const double);
    inline double grad4xbiastop(const double, const double, const double, const double);
    */
};
#endif
