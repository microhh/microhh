#ifndef MPIINTERFACE
#define MPIINTERFACE

#include <mpi.h>
#include "input.h"
#include "grid.h"

class cmpi
{
  public:
    cmpi(cgrid *);
    ~cmpi();

    int readinifile(cinput *);
    int init();

    int boundary_cyclic(double *);
    int transposezx(double *, double *);
    int transposexz(double *, double *);

    int nprocs;
    int npx;
    int npy;
    int mpiid;
    int mpicoordx;
    int mpicoordy;

    int nnorth;
    int nsouth;
    int neast;
    int nwest;

  private:
    cgrid *grid;

    bool initialized;

    MPI_Comm commxy;
    MPI_Comm commx;
    MPI_Comm commy;

    MPI_Datatype eastwestedge;
    MPI_Datatype northsouthedge;
    MPI_Datatype transposez;
    MPI_Datatype transposex;

    MPI_Request *reqsx;
    MPI_Request *reqsy;
};
#endif
