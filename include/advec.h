#ifndef ADVEC
#define ADVEC

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cadvec
{
  public:
    cadvec(cgrid *, cfields *, cmpi *);
    ~cadvec();

    int exec();
    double getcfl(double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double calccfl_2nd(double *, double *, double *, double *, double);
    int advecu_2nd(double *, double *, double *, double *, double *);
    int advecv_2nd(double *, double *, double *, double *, double *);
    int advecw_2nd(double *, double *, double *, double *, double *);
    int advecs_2nd(double *, double *, double *, double *, double *, double *);
    inline double interp2(const double, const double);

    double calccfl_4th(double *, double *, double *, double *, double);
    int advecu_4th(double *, double *, double *, double *, double *);
    int advecv_4th(double *, double *, double *, double *, double *);
    int advecw_4th(double *, double *, double *, double *, double *);
    int advecs_4th(double *, double *, double *, double *, double *, double *);
    inline double interp4(const double, const double, const double, const double);
    inline double interp4bot(const double, const double, const double, const double);
    inline double interp4top(const double, const double, const double, const double);
};
#endif
