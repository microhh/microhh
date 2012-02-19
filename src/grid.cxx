#include "grid.h"

// build the grid
grid::grid()
{
  int k;

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

  dz    = new double[kmax];
  dzh   = new double[kmax];

  dx    = xsize / itot;
  dy    = ysize / jtot;

  for(k=0; k<kmax; k++)
  {
    dz [k] = zsize / ktot;
    dzh[k] = zsize / ktot;
  }
}

