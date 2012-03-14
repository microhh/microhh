#ifndef PRES
#define PRES

#include "grid.h"
#include "fields.h"

class cpres
{
  public:
    cpres(cgrid *, cfields *);
    ~cpres();

    int exec(double);

  private:
    cgrid *grid;
    cfields *fields;

    int pres_2nd_in(double *, 
                    double *, double *, double *,
                    double *, double *, double *,
                    double *, double);
    int pres_2nd_out(double *, double *, double *,
                     double *, double *);
};
#endif
