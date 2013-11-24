#ifndef THERMO_MOIST
#define THERMO_MOIST

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "thermo.h"
#include <cmath>

#define rd 287.04
#define rv 461.5
#define ep rd/rv
#define cp 1005
#define lv 2.5e6
#define rhow    1.e3
#define tmelt   273.15
#define p0 1.e5
#define grav 9.81

class cthermo_moist : public cthermo
{
  public:
    cthermo_moist(cgrid *, cfields *, cmpi *);
    ~cthermo_moist();
    int readinifile(cinput *);
    int init(cinput *);
    int create();
    int exec();
    int getbuoyancy();
    int getbuoyancyh();
    int getsat();
    int getsath();

  private:
    double ps;
    double thvs;
    double rhos;
    double *pmn;

    bool allocated;
    
    int buoyancy_2nd(double *, double *, double *, double *);
    int buoyancy_4th(double *, double *, double *, double *);
    int calcqlfield(double *, double *, double *, double *);
    inline double calcql(const double, const double, const double);
    inline double bu(const double p, const double s, const double qt, const double ql);
    inline double exner(const double);
    inline double rslf(const double, const double);
    inline double esl(const double);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
    inline double interp4biasbot(const double, const double, const double, const double);
    inline double interp4biastop(const double, const double, const double, const double);

};
#endif
