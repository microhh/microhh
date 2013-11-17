#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "thermo_moist.h"
#include "defines.h"

cthermo_moist::cthermo_moist(cgrid *gridin, cfields *fieldsin, cmpi *mpiin) : cthermo(gridin, fieldsin, mpiin)
{
  allocated = false;
}

cthermo_moist::~cthermo_moist()
{
  // std::printf("Destroying instance of object buoyancy\n");
  if (allocated)
  {
    delete[] pmn;
  }
}

int cthermo_moist::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&ps, "thermo", "ps", "");

  nerror += fields->initpfld("s");
  nerror += inputin->getItem(&fields->sp["s"]->visc, "fields", "svisc", "s");
  nerror += fields->initpfld("qt");
  nerror += inputin->getItem(&fields->sp["qt"]->visc, "fields", "svisc", "qt");
  nerror += fields->initdfld("ql");

  return (nerror > 0);

}

int cthermo_moist::create()
{
  int nerror = 0;
  double qtsurf, ssurf;
  
  // Create hydrostatic profile
  nerror += grid->calcmean(&ssurf, fields->s["s"]->databot,1);
  nerror += grid->calcmean(&qtsurf, fields->s["qt"]->databot,1);

  thvs = ssurf * (1. - (1. - rv/rd)*qtsurf);
  double tvs  = exner(ps) * thvs;
  rhos = ps / (rv * tvs);

  pmn = new double[grid->kcells];
  for(int k=0; k<grid->kcells; k++)
  {
    pmn[k] = ps - rhos * grav * grid->z[k];
  }

  allocated = true;
  return nerror;
}
int cthermo_moist::exec()
{

  // add mean pressure to pressure fluctuations into tmp array
  for(int k=0; k<grid->kcells; k++)
  {
    for(int n=0; n<grid->icells*grid->jcells; n++)
    {
      fields->s["tmp1"]->data[n] = (fields->s["p"]->data[n] + pmn[k]) / rhos;
    }
  }


  // extend later for gravity vector not normal to surface
  if(grid->swspatialorder == "2")
  {
    buoyancy_2nd(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, fields->s["tmp1"]->data);
  }
  else if(grid->swspatialorder == "4")
  {
    buoyancy_4th(fields->wt->data, fields->s["s"]->data, fields->s["qt"]->data, fields->s["tmp1"]->data);
  }

  return 0;
}

int cthermo_moist::buoyancy_2nd(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj,kk;
  double sh, qth, ph, ql;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // CvH check the usage of the gravity term here, in case of scaled DNS we use one. But thermal expansion coeff??
  for(int k=grid->kstart+1; k<grid->kend; k++)
  {
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        sh  = interp2(s[ijk-kk], s[ijk]);
        qth = interp2(qt[ijk-kk], qt[ijk]);
        ph  = interp2(p[ijk-kk], p[ijk]);
        ql  = calcql(sh, qth, ph);
        wt[ijk] += bu(sh, qth, ql);
      }
  }
  return 0;
}

int cthermo_moist::buoyancy_4th(double * restrict wt, double * restrict s, double * restrict qt, double * restrict p)
{
  int ijk,jj;
  int kk1,kk2;
  double sh, qth, ph, ql;

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
        ql  = calcql(sh, qth, ph);
        wt[ijk] += bu(sh, qth, ql);
      }

  return 0;
}

inline double cthermo_moist::bu(const double s, const double qt, const double ql)
{
  return grav * (s * (1. - (1. - rv/rd)*qt - rv/rd*ql) - thvs) / thvs;
}

inline double cthermo_moist::calcql(const double s, const double qt, const double p)
{
  int niter = 0, nitermax = 5;
  double sabs, sguess = 1.e9, t, ql, dtldt;
  sabs = s * exner(p);
  t = sabs;
  while (fabs(sguess-sabs)/sabs > 1e-5 && niter < nitermax)
  {
    ++niter;
    ql = std::max(0.,qt - rslf(p, t));
    sguess = t*exp(-lv*ql/(cp*t));
    dtldt = sabs/t*(1. - lv * ql / (pow(cp,2.)* t));
    t += (sguess-sabs) / dtldt;
  }
  return ql;
}

inline double cthermo_moist::exner(const double p)
{
  return pow((p/p0),(rd/cp));
}
inline double cthermo_moist::rslf(const double p, const double t)
{
  return ep*esl(t)/(p-esl(t));
}

inline double cthermo_moist::esl(const double t)
{
  const double c0=0.6105851e+03;
  const double c1=0.4440316e+02;
  const double c2=0.1430341e+01;
  const double c3=0.2641412e-01;
  const double c4=0.2995057e-03;
  const double c5=0.2031998e-05;
  const double c6=0.6936113e-08;
  const double c7=0.2564861e-11;
  const double c8=-.3704404e-13;
  const double x=std::max(-80.,t-tmelt);

  return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
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

