#ifndef FIELDS
#define FIELDS

#include "grid.h"

class cfield3d
{
  public:
    // functions
    cfield3d(cgrid *, double *);
    ~cfield3d();
    int index(int, int, int);
    int boundary_bottop(int);
    int boundary_cyclic();

    // variables
    double *data;

  private:
    cgrid *grid;
};

class cfields
{
  public:
    // functions
    cfields(cgrid *);
    ~cfields();
    int initfields();
    int createfields();

    int resettend();
    int boundary_bottop();

    // variables
    cfield3d *u;
    cfield3d *v;
    cfield3d *w;

    cfield3d *ut;
    cfield3d *vt;
    cfield3d *wt;

  private:
    // variables
    cgrid *grid;
    double *flow;
    double *flowt;
};
#endif

