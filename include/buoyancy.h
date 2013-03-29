#ifndef BUOYANCY
#define BUOYANCY

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cbuoyancy
{
  public:
    cbuoyancy(cgrid *, cfields *, cmpi *);
    ~cbuoyancy();
    int readinifile(cinput *);
    int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    std::string swbuoyancy;
    double gravitybeta; // gravity multiplied with thermal expansion coefficient

    int buoyancy_2nd(double *, double *);
    int buoyancy_4th(double *, double *);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
    inline double interp4biasbot(const double, const double, const double, const double);
    inline double interp4biastop(const double, const double, const double, const double);
};
#endif
