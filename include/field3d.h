#ifndef FIELD3D
#define FIELD3D

#include <string>
#include "grid.h"

class cfield3d
{
  public:
    // functions
    cfield3d(cgrid *, std::string);
    ~cfield3d();
    int init();
    int boundary_bottop(int);
    int boundary_cyclic();
    int save(int, int);
    int load(int, int);

    // variables
    double *data;
    std::string name;

  private:
    cgrid *grid;
    bool allocated;
};
#endif

