#ifndef ADVEC_G4
#define ADVEC_G4

#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "mpiinterface.h"

class cadvec_g4 : public cadvec
{
  public:
    cadvec_g4(cgrid *, cfields *, cmpi *);
    ~cadvec_g4();

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
};
#endif
