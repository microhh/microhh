#include <iostream>
#include <cstdio>
#include <cmath>
#include "grid.h"

// build the grid
grid::grid()
{
  std::cout << "Creating grid" << std::endl;
  xsize = 6.28;
  ysize = 3.14;
  zsize = 2.;

  itot  = 32;
  jtot  = 32;
  ktot  = 32;

  igc   = 1;
  jgc   = 1;
  kgc   = 1;

  imax  = itot;
  jmax  = jtot;
  kmax  = ktot;

  // convenience variables
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

  dx    = xsize / itot;
  dy    = ysize / jtot;

  // create non-equidistant grid
  double alpha = 0.967;
  double eta;

  int k;
  for(k=kstart; k<kend; k++)
  {
    eta  = -1. + 2.*((k-kstart+1) - 0.5) / kmax;
    z[k] = zsize / (2.*alpha) * std::tanh(eta*0.5*(std::log(1.+alpha) - std::log(1.-alpha))) + 0.5*zsize;
    //z[k] = zsize / (2*kmax) + zsize / kmax * (k-kstart);
  }

  for(k=0; k<kgc; k++)
  {
    z[kstart-k-1] = -1. * z[kstart+k];
    z[kend  +k  ] = -1. * z[kend-1-k] + 2.*zsize;
  }

  //do k = 2-kgc, kmax+kgc
  for(k=kstart-kgc+1; k<kend+kgc; k++)
  {
    dzh[k] = z[k] - z[k-1];
    zh [k] = 0.5*(z[k] + z[k-1]);
  }

  for(k=kstart; k<kend; k++)
    dz[k]  = 0.5*(z[k]-z[k-1]) + 0.5*(z[k+1]-z[k]);

  for(k=0; k<kgc; k++)
  {
    dz[kstart-k-1] = dz[kstart+k];
    dz[kend+k]     = dz[kend-k-1];
  }

  //for(int k=kstart-kgc; k<kend+kgc; k++)
  //  std::printf("%4d %9.6f %9.6f %9.6f %9.6f \n", k-kstart+1, z[k], zh[k], dz[k], dzh[k]);
}

grid::~grid()
{
  std::cout << "Deleting grid" << std::endl;
  delete[] dz;
  delete[] dzh;
}

