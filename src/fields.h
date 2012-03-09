#ifndef FIELDS
#define FIELDS

#include "grid.h"

class cfield
{
  public:
    // functions
    cfield(cgrid *, double *);
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
    int resettend();

    // variables
    cfield *u;
    cfield *v;
    cfield *w;

    cfield *ut;
    cfield *vt;
    cfield *wt;
  private:
    // variables
    cgrid *grid;
    double *flow;
    double *flowt;
};
#endif
