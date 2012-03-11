#ifndef ADVEC
#define ADVEC

#include "grid.h"
#include "fields.h"

class cadvec
{
  public:
    cadvec(cgrid *, cfields *);
    ~cadvec();

    int exec(double);

  private:
    cgrid *grid;
    cfields *fields;

    double courant(double *, double *, double *, double *, double);
    int advecu_2nd(double *, double *, double *, double *, double *);
    int advecv_2nd(double *, double *, double *, double *, double *);
    int advecw_2nd(double *, double *, double *, double *, double *);
    inline double interp2(const double, const double);
};
#endif
