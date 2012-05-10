#ifndef GRID
#define GRID

#include <mpi.h>
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
    int create();
    int calculate();
    int save(int);
    int load(int);

    // variables
    int itot;
    int jtot;
    int ktot;

    int imax;
    int jmax;
    int kmax;

    int iblock;
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
    
    double *x;
    double *y;
    double *z;
    double *xh;
    double *yh;
    double *zh;

    // MPI functions
    int initmpi();
    int boundary_cyclic(double *);
    int transposezx(double *, double *);
    int transposexz(double *, double *);
    int transposexy(double *, double *);
    int transposeyx(double *, double *);
    int transposeyz(double *, double *);

    int getmax(double *);
    int getsum(double *);

    int writefield3d(double *, char *);
    int readfield3d (double *, char *);

  private:
    cmpi *mpi;
    bool allocated;

    // MPI Datatypes
    MPI_Datatype eastwestedge;
    MPI_Datatype northsouthedge;
    MPI_Datatype transposez;
    MPI_Datatype transposex;
    MPI_Datatype transposex2;
    MPI_Datatype transposey;
    MPI_Datatype subarray;
};
#endif
