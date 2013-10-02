#ifndef FORCE
#define FORCE

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cforce
{
  public:
    cforce(cgrid *, cfields *, cmpi *);
    ~cforce();
    int readinifile(cinput *);
    int init();
    int create(cinput *);
    int exec(double);

    int save();
    int load();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    bool allocated;

    std::string swforce;

    double uflow;
    double fc;

    double *ug;
    double *vg;

    int flux(double *, double *, double *, double);
    int coriolis_2nd(double *, double *, double *, double *, double *, double *);
    int coriolis_4th(double *, double *, double *, double *, double *, double *);
};
#endif
