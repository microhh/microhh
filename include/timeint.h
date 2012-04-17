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
    double getsubdt(double);

    // variables
    int substep;

  private:
    cgrid *grid;
    cfields *fields;

    int rk3(double *, double *, double);
    int rk4(double *, double *, double);
    double rk3subdt(double);
    double rk4subdt(double);
};
#endif
