#ifndef GRID
#define GRID

#ifdef PARALLEL
#include <mpi.h>
#endif
#include <fftw3.h>
#include "input.h"
#include "mpiinterface.h"

class cgrid
{
  public:
    // functions
    cgrid(cmpi *);
    ~cgrid();
    int readinifile(cinput *);
    int init();
    int create(cinput *);
    int calculate();
    int save();
    int load();

    // variables
    int itot;
    int jtot;
    int ktot;

    int imax;
    int jmax;
    int kmax;

    int iblock;
    int jblock;
    int kblock;

    int igc;
    int jgc;
    int kgc;

    int icells;
    int jcells;
    int kcells;
    int ncells;
    int istart;
    int jstart;
    int kstart;
    int iend;
    int jend;
    int kend;

    double xsize;
    double ysize;
    double zsize;

    double dx;
    double dy;
    double *dz;
    double *dzh;
    double *dzi;
    double *dzhi;
    double *dzi4;
    double *dzhi4;
    
    double *x;
    double *y;
    double *z;
    double *xh;
    double *yh;
    double *zh;

    std::string swspatialorder;

    // MPI functions
    int initmpi();
    int exitmpi();
    int boundary_cyclic  (double *);
    int boundary_cyclic2d(double *);
    int transposezx(double *, double *);
    int transposexz(double *, double *);
    int transposexy(double *, double *);
    int transposeyx(double *, double *);
    int transposeyz(double *, double *);
    int transposezy(double *, double *);

    int getmax (double *);
    int getsum (double *);
    int getprof(double *, int);

    int savefield3d(double *, double *, double *, char *);
    int loadfield3d(double *, double *, double *, char *);

    int savexzslice(double *, double *, int, char *);
    int savexyslice(double *, double *, int, char *);

    // variables for the fast fourier transforms
    double *fftini, *fftouti;
    double *fftinj, *fftoutj;
    fftw_plan iplanf, iplanb;
    fftw_plan jplanf, jplanb;

    int fftforward (double *, double *, double *, double *, double *, double *);
    int fftbackward(double *, double *, double *, double *, double *, double *);

    // interpolation functions
    int interpolatex_2nd(double *, double *, int);
    int interpolatey_2nd(double *, double *, int);
    int interpolatex_4th(double *, double *, int);
    int interpolatey_4th(double *, double *, int);
    // int interpolatez_4th(double *, double *, int);

  private:
    cmpi *mpi;
    bool allocated;
    bool mpitypes;
    bool fftwplan;

#ifdef PARALLEL
    // MPI Datatypes
    MPI_Datatype eastwestedge;
    MPI_Datatype northsouthedge;
    MPI_Datatype eastwestedge2d;
    MPI_Datatype northsouthedge2d;
    MPI_Datatype transposez;
    MPI_Datatype transposez2;
    MPI_Datatype transposex;
    MPI_Datatype transposex2;
    MPI_Datatype transposey;
    MPI_Datatype transposey2;
    MPI_Datatype subi;
    MPI_Datatype subj;
    MPI_Datatype subarray;
    MPI_Datatype subxzslice;
    MPI_Datatype subxyslice;
    double *profl;
#endif
};
#endif
