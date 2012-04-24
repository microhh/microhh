#ifndef MPIINTERFACE
#define MPIINTERFACE

#include <mpi.h>
// #include "grid.h"
// #include "fields.h"
#include "input.h"

class cmpi
{
  public:
    // cmpi(cgrid *, cfields *);
    cmpi();
    ~cmpi();

    int readinifile(cinput *);
    int init();

  private:
    int npx;
    int npy;
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
