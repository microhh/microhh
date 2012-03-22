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
    bool insubstep();
    double subdt(double);

    // variables
    int substep;

  private:
    cgrid *grid;
    cfields *fields;

    int rk3(double *, double *, double *, double *, double);
    int rk4(double *, double *, double *, double *, double);
};
#endif
