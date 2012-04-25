#ifndef MPIINTERFACE
#define MPIINTERFACE

#include <mpi.h>
#include "input.h"

class cmpi
{
  public:
    cmpi();
    ~cmpi();

    int readinifile(cinput *);
    int init(cgrid *);

    int boundary_cyclic(double *, cgrid *);

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
    MPI_Comm commxy;
    MPI_Comm commx;
    MPI_Comm commy;

    MPI_Datatype eastwestedge;
    MPI_Datatype northsouthedge;
};
#endif
