#ifndef MPIINTERFACE
#define MPIINTERFACE

#ifdef PARALLEL
#include <mpi.h>
#endif
#include <string>
#include "input.h"

class cmpi
{
  public:
    cmpi();
    ~cmpi();

    int startup(int, char**);
    int readinifile(cinput *);
    int init();

    double gettime();
    int waitall();

    // overload the broadcast function
    int broadcast(char *, int);
    int broadcast(int *, int);
    int broadcast(double *, int);
    int broadcast(unsigned long *, int);

    std::string mode;
    std::string simname;
    std::string sw_model;

    int nprocs;
    int npx;
    int npy;
    int mpiid;
    int mpicoordx;
    int mpicoordy;

#ifdef PARALLEL
    int nnorth;
    int nsouth;
    int neast;
    int nwest;

    MPI_Comm commxy;
    MPI_Comm commx;
    MPI_Comm commy;

    MPI_Request *reqs;
    int reqsn;
#endif

  private:
    bool initialized;
    bool allocated;

#ifdef PARALLEL
    int checkerror(int);
#endif
};
#endif
