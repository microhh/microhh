#ifndef TIMEINT
#define TIMEINT

#include "grid.h"
#include "fields.h"

class ctimeint
{
  public:
    // functions
    ctimeint(cgrid *, cfields *);
    ~ctimeint();

    int exec(double);
    int rk3(double *, double *, double);
    bool insubstep();

    // variables
    int substep;

  private:
    cgrid *grid;
    cfields *fields;
};
#endif
