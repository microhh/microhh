#include <cstdio>
#include <cmath>
#include "grid.h"

// build the grid
cgrid::cgrid()
{
  std::printf("Creating instance of object grid\n");
  xsize = 6.28;
  ysize = 3.14;
  zsize = 2.;

  itot  = 32;
  jtot  = 32;
  ktot  = 32;

  igc   = 1;
  jgc   = 1;
  kgc   = 1;
}

cgrid::~cgrid()
{
  delete[] dz;
  delete[] dzh;
  delete[] dzi;
  delete[] dzhi;
  std::printf("Destroying instance of object grid\n");
}

int cgrid::initgrid()
{
  imax  = itot;
  jmax  = jtot;
  kmax  = ktot;

  icells = (imax+2*igc);
  jcells = (jmax+2*jgc);
  kcells = (kmax+2*kgc);
  ncells = (imax+2*igc)*(jmax+2*jgc)*(kmax+2*kgc);

  istart = igc;
  jstart = jgc;
  kstart = kgc;

  iend   = imax + igc;
  jend   = jmax + jgc;
  kend   = kmax + kgc;

  z     = new double[kmax+2*kgc];
  zh    = new double[kmax+2*kgc];
  dz    = new double[kmax+2*kgc];
  dzh   = new double[kmax+2*kgc];
  dzi   = new double[kmax+2*kgc];
  dzhi  = new double[kmax+2*kgc];

  dx    = xsize / itot;
  dy    = ysize / jtot;

  return 0;
}

int cgrid::creategrid()
{
  // create non-equidistant grid
  double alpha = 0.967;
  double eta;
  int k;
  // heights are set according to Moser180 case
  for(k=kstart; k<kend; k++)
  {
    eta  = -1. + 2.*((k-kstart+1) - 0.5) / kmax;
    z[k] = zsize / (2.*alpha) * std::tanh(eta*0.5*(std::log(1.+alpha) - std::log(1.-alpha))) + 0.5*zsize;
    //z[k] = zsize / (2*kmax) + zsize / kmax * (k-kstart);
  }
  // end Moser180 setup

  // calculate the height of the ghost cells
  for(k=0; k<kgc; k++)
  {
    z[kstart-k-1] = -1. * z[kstart+k];
    z[kend  +k  ] = -1. * z[kend-1-k] + 2.*zsize;
  }

  // assume the flux levels are exactly in between the cells
  // compute the flux levels and the distance between them
  for(k=1; k<kcells; k++)
  {
    zh  [k] = 0.5*(z[k] + z[k-1]);
    dzh [k] = z[k] - z[k-1];
    dzhi[k] = 1./dzh[k];
  }

  // compute the heigth of the grid cells
  for(k=kstart; k<kend; k++)
  {
    dz [k] = 0.5*(z[k]-z[k-1]) + 0.5*(z[k+1]-z[k]);
    dzi[k] = 1./dz[k];
  }

  // compute the height of the ghost cells
  for(k=0; k<kgc; k++)
  {
    dz[kstart-k-1] = dz[kstart+k];
    dz[kend+k]     = dz[kend-k-1];
  }

  return 0;
}
