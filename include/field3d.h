#ifndef FIELD3D
#define FIELD3D

#include <string>
#include "grid.h"
#include "mpiinterface.h"

class cfield3d
{
  public:
    // functions
    cfield3d(cgrid *, cmpi *, std::string);
    ~cfield3d();

    int init();
    // int boundary_bottop(int);
    // int boundary_cyclic();
    int save(int, double *);
    int load(int, double *);

    // variables
    double *data;
    std::string name;

  private:
    cgrid *grid;
    cmpi  *mpi;
    bool allocated;
};
#endif

