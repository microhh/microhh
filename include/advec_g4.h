#ifndef ADVEC_G4
#define ADVEC_G4

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "advec.h"

class cadvec_g4 : public cadvec
{
  public:
    cadvec_g4(cgrid *, cfields *, cmpi *);
    ~cadvec_g4();

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
};
#endif
