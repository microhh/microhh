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
    int getql(cfield3d *, cfield3d *);

    // functions to retrieve buoyancy properties, to be called from other classes
    int getbuoyancysurf   (cfield3d *);
    int getbuoyancyfluxbot(cfield3d *);
    int getbuoyancy(cfield3d *, cfield3d *);
    // int getbuoyancyh();
    // int getsat();
    // int getsath();

  private:
    double ps;
    double thvs;
    double rhos;
    double *pmn;

    bool allocated;
    
    int calcbuoyancytend_2nd(double *, double *, double *, double *);
    int calcbuoyancytend_4th(double *, double *, double *, double *);

    int calcbuoyancy(double *, double *, double *, double *);

    int calcpres(double *, double *, double *);
    int calcqlfield(double *, double *, double *, double *);
    int calcbuoyancybot(double *, double *,
                        double *, double *,
                        double *, double *);
    int calcbuoyancyfluxbot(double *, double *, double *, double *, double *);

    inline double calcql(const double, const double, const double);
    inline double bu(const double, const double, const double, const double);
    inline double bunoql(const double, const double, const double);
    inline double bufluxnoql(const double, const double, const double, const double, const double);
    inline double exner(const double);
    inline double rslf(const double, const double);
    inline double esl(const double);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
};
#endif
