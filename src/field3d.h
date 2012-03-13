#ifndef FIELD3D
#define FIELD3D

#include "grid.h"

class cfield3d
{
  public:
    // functions
    cfield3d(cgrid *, double *);
    ~cfield3d();
    int boundary_bottop(int);
    int boundary_cyclic();

    // variables
    double *data;

  private:
    cgrid *grid;
};
#endif

