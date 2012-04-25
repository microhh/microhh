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
    int checkLayout();
    int checkBoundary();
    int checkTranspose();

  private:
    cgrid   *grid;
    cmpi    *mpi;

    cfield3d *s;
    cfield3d *temp1;
    cfield3d *temp2;
};
#endif
