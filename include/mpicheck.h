#ifndef MPICHECK
#define MPICHECK

#include "grid.h"
#include "field3d.h"
#include "mpiinterface.h"

class cmpicheck
{
  public:
    cmpicheck(cgrid *, cmpi *);
    ~cmpicheck();
    
    int create();
    int boundary();

    int showLayout();
    int showLine();

  private:
    cgrid   *grid;
    cmpi    *mpi;

    cfield3d *s;
};
#endif
