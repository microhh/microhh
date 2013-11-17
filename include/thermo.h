#ifndef THERMO
#define THERMO

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cthermo
{
  public:
    cthermo(cgrid *, cfields *, cmpi *);
    virtual ~cthermo();
    virtual int readinifile(cinput *);
    virtual int exec();

  protected:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

//     std::string swthermo;
    double gravitybeta; // gravity multiplied with thermal expansion coefficient

  private:
    int buoyancy_2nd(double *, double *);
    int buoyancy_4th(double *, double *);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
    inline double interp4biasbot(const double, const double, const double, const double);
    inline double interp4biastop(const double, const double, const double, const double);
};
#endif
