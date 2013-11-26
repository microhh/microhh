#ifndef THERMO_DRY
#define THERMO_DRY

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "thermo.h"

/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */

class cthermo_dry : public cthermo
{
  public:
    cthermo_dry(cgrid *, cfields *, cmpi *); ///< Constructor of the dry thermodynamics class.
    ~cthermo_dry();                          ///< Destructor of the dry thermodynamics class.
    int readinifile(cinput *);               ///< Processing data of the input file.
    int exec();                              ///< Add the tendencies belonging to the buoyancy.

    int getbuoyancy(cfield3d *, cfield3d *); ///< Compute the buoyancy for usage in another routine.
    int getbuoyancyfluxbot(cfield3d *);      ///< Compute the bottom buoyancy flux for usage in another routine.

  private:
    int calcbuoyancy(double *, double *);         ///< Calculation of the buoyancy.
    int calcbuoyancyfluxbot(double *, double *);  ///< Calculation of the buoyancy flux at the bottom.
    int calcbuoyancytend_2nd(double *, double *); ///< Calculation of the buoyancy tendency with 2nd order accuracy.
    int calcbuoyancytend_4th(double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.

    inline double interp2(const double, const double); ///< 2nd order interpolation function.
    inline double interp4(const double, const double, 
                          const double, const double); ///< 4th order interpolation function.

    double gravitybeta; ///< Gravity acceleration multiplied with thermal expansion coefficient.
};
#endif
