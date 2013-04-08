#ifndef ADVEC_G2I4
#define ADVEC_G2I4

#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "mpiinterface.h"

class cadvec_g2i4 : public cadvec
{
  public:
    cadvec_g2i4(cgrid *, cfields *, cmpi *);
    ~cadvec_g2i4();

    unsigned long gettimelim(long unsigned int idt, double ifactor);
    double getcfl(double);
    int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double calccfl(double *, double *, double *, double *, double);
    int advecu(double *, double *, double *, double *, double *);
    int advecv(double *, double *, double *, double *, double *);
    int advecw(double *, double *, double *, double *, double *);
    int advecs(double *, double *, double *, double *, double *, double *);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
    inline double interp4bot(const double, const double, const double, const double);
    inline double interp4top(const double, const double, const double, const double);
};
#endif
