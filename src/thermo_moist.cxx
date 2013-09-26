#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "defines.h"

cthermo_moist::cthermo_moist(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object buoyancy\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cthermo_moist::~cthermo_moist()
{
  // std::printf("Destroying instance of object buoyancy\n");
}

int cthermo_moist::readinifile(cinput *inputin)
{

}

int cthermo_moist::init(cinput *inputin)
{
  int nerror=0;

  nerror += fields->initpfld("s");
  nerror += inputin->getItem(&fields->sp["s"]->visc, "fields", "svisc", "s");
  nerror += fields->initpfld("qt");
  nerror += inputin->getItem(&fields->sp["qt"]->visc, "fields", "svisc", "qt");
  nerror += fields->initdfld("ql");

}

int cthermo_moist::exec()
{
  if(swbuoyancy == "0")
    return 0;

  // extend later for gravity vector not normal to surface
  if(swbuoyancy == "2")
    buoyancy_2nd(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, fields->s["p"]->data);
  else if(swbuoyancy == "4")
    buoyancy_4th(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, fields->s["p"]->data);

  return 0;
}

int cthermo_moist::buoyancy_2nd(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj,kk;
  double sh, qth, ph;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // CvH check the usage of the gravity term here, in case of scaled DNS we use one. But thermal expansion coeff??
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        sh  = interp2(s[ijk-kk], s[ijk]);
        qth = interp2(qt[ijk-kk], qt[ijk]);
        ph  = interp2(p[ijk-kk], p[ijk]);

        wt[ijk] += bu(sh, qth,calcql(sh,qth,ph));
      }

  return 0;
}

int cthermo_moist::buoyancy_4th(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj;
  int kk1,kk2;
  double sh, qth, ph;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk1;
        sh  = interp4(s[ijk-kk2] , s[ijk-kk1] , s[ijk] , s[ijk+kk1]);
        qth = interp4(qt[ijk-kk2], qt[ijk-kk1], qt[ijk], qt[ijk+kk1]);
        ph  = interp4(p[ijk-kk2] , p[ijk-kk1] , p[ijk] , p[ijk+kk1]);

        wt[ijk] += bu(sh, qth, calcql(sh,qth,ph));
      }

  return 0;
}

inline double cthermo_moist::bu(const double s, const double qt, const double ql)
{
  return gravitybeta * s * (1. - (1. - rv/rd)*qt - rv/rd*ql);
}

inline double cthermo_moist::calcql(const double s, const double qt, const double p)
{
  return 0.;
}
inline double cthermo_moist::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cthermo_moist::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

inline double cthermo_moist::interp4biasbot(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*a + (15./16.)*b - (5./16.)*c + (1./16)*d);
}

inline double cthermo_moist::interp4biastop(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*d + (15./16.)*c - (5./16.)*b + (1./16)*a);
}

