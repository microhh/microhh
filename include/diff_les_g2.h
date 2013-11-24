#ifndef DIFF_LES_G2
#define DIFF_LES_G2

#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "mpiinterface.h"
#include "boundary_surface.h"
#include "thermo.h"

class cdiff_les_g2 : public cdiff
{
  public:
    cdiff_les_g2(cgrid *, cfields *, cmpi *);
    ~cdiff_les_g2();

    int readinifile(cinput *);
    int exec();
    int execvisc();

    unsigned long gettimelim(unsigned long, double);
    double getdn(double);

    int setdepends(cthermo *, cboundary_surface *);

  private:
    int strain2(double *,
                double *, double *, double *,
                double *, double *,
                double *, double *,
                double *, double *, double *);
    int evisc(double *,
              double *, double *, double *, double *,
              double *, double *, double *,
              double *, double *,
              double *, double *, double *,
              double);
    int evisc_neutral(double *,
                      double *, double *, double *,
                      double *, double *,
                      double *, double *);
    int diffu(double *, double *, double *, double *, double *, double *, double *, double *, double *);
    int diffv(double *, double *, double *, double *, double *, double *, double *, double *, double *);
    int diffw(double *, double *, double *, double *, double *, double *, double *);
    int diffc(double *, double *, double *, double *, double *, double *, double *, double);

    double calcdnmul(double *, double *, double);

    inline double phim(double);
    inline double phih(double);

    double cs;

    cthermo           *thermo;
    cboundary_surface *boundary;
};
#endif
