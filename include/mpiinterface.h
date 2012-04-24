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
    int init();
    int npx;
    int npy;

  private:
    int nprocs;
    int mpiid;
    int mpicoordx;
    int mpicoordy;

    int nnorth;
    int nsouth;
    int neast;
    int nwest;

    MPI_Comm commxy;
    MPI_Comm commx;
    MPI_Comm commy;
};
#endif
