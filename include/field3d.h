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
    int save(int, double *, double *);
    int load(int, double *, double *);
    int checkfornan();

    // variables
    double *data;
    double *databot;
    double *datatop;
    double *datagradbot;
    double *datagradtop;
    double *datafluxbot;
    double *datafluxtop;
    std::string name;
    double visc;

  private:
    cgrid *grid;
    cmpi  *mpi;
    bool allocated;
};
#endif

