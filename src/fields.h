#ifndef FIELDS
#define FIELDS

#include "grid.h"

class cfields
{
  public:
    // functions
    cfields(cgrid*);
    ~cfields();
    int resettend();

    // variables
    double *u;
    double *v;
    double *w;

    double *ut;
    double *vt;
    double *wt;
  private:
    // variables
    cgrid *grid;
    double *flow;
    double *flowt;
};
#endif
