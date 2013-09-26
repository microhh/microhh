#ifndef THERMO_MOIST
#define THERMO_MOIST

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cthermo_moist
{
  public:
    cthermo_moist(cgrid *, cfields *, cmpi *);
    ~cthermo_moist();
    int readinifile(cinput *);
    int init(cinput *);
    int exec();
    int getbuoyancy();
    int getbuoyancyh();
    int getsat();
    int getsath();

#define rd 287.04
#define rv 461.5
#define rm 461.5
#define ep rd/rm
#define cp 1005
#define latheat 2.5e6
#define rhow    1.e3
    
  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    std::string swbuoyancy;
    double gravitybeta; // gravity multiplied with thermal expansion coefficient

    int buoyancy_2nd(double *, double *, double *, double *);
    int buoyancy_4th(double *, double *, double *, double *);
    inline double calcql(const double, const double, const double);
    inline double bu(const double, const double, const double);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
    inline double interp4biasbot(const double, const double, const double, const double);
    inline double interp4biastop(const double, const double, const double, const double);

};
#endif