#ifndef FIELDS
#define FIELDS

#include "grid.h"

class fields
{
  public:
    // functions
    fields(grid*);
    ~fields();

    // variables
    double *u;
    double *v;
    double *w;

    double *ut;
    double *vt;
    double *wt;
  private:
    // variables
    double *flow;
};
#endif
