#ifndef THERMO_DRY
#define THERMO_DRY

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "thermo.h"

class cthermo_dry : public cthermo
{
  public:
    cthermo_dry(cgrid *, cfields *, cmpi *);
    ~cthermo_dry();
    int readinifile(cinput *);
    int exec();

  private:
    int buoyancy_2nd(double *, double *);
    int buoyancy_4th(double *, double *);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
    inline double interp4biasbot(const double, const double, const double, const double);
    inline double interp4biastop(const double, const double, const double, const double);
};
#endif
