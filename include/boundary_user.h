#ifndef BOUNDARY_USER
#define BOUNDARY_USER

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "boundary.h"

class cboundary_user : public cboundary
{
  public:
    cboundary_user(cgrid *, cfields *, cmpi *);

    int readinifile(cinput *);
    int setvalues();

  private:
    int setbc_patch(double *, double, double, double);

    // patch type
    int    patch_dim;
    double patch_xh;
    double patch_xr;
    double patch_xi;
    double patch_facr;
    double patch_facl;
};
#endif
