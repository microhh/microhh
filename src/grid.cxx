#include <cstdio>
#include <cmath>
#include "grid.h"
#include "input.h"
#include "defines.h"

// build the grid
cgrid::cgrid(cmpi *mpiin)
{
  // std::printf("Creating instance of object grid\n");
  mpi = mpiin;

  allocated = false;
  mpitypes  = false;
}

cgrid::~cgrid()
{
  if(allocated)
  { 
    delete[] x;
    delete[] xh;
    delete[] y;
    delete[] yh;
    delete[] z;
    delete[] zh;
    delete[] dz;
    delete[] dzh;
    delete[] dzi;
    delete[] dzhi;
    delete[] dzi4;
    delete[] dzhi4;
  }

  exitmpi();
}

int cgrid::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&xsize, "grid", "xsize");
  n += inputin->getItem(&ysize, "grid", "ysize");
  n += inputin->getItem(&zsize, "grid", "zsize");

  n += inputin->getItem(&itot, "grid", "itot");
  n += inputin->getItem(&jtot, "grid", "jtot");
  n += inputin->getItem(&ktot, "grid", "ktot");

  if(n > 0)
    return 1;
  
  igc = 3;
  jgc = 3;
  kgc = 3;

  return 0;
}

int cgrid::init()
{
  // check whether the grid fits the processor configuration
  if(itot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npx = %d\n", itot, mpi->npx);
    return 1;
  }
  if(itot % mpi->npy != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npy = %d\n", itot, mpi->npy);
    return 1;
  }
  // check this one only when npy > 1, since the transpose in that direction only happens then
  //if(jtot % mpi->npx != 0 && mpi->npy > 1)
  if(jtot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npx = %d\n", jtot, mpi->npx);
    return 1;
  }
  if(jtot % mpi->npy != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npy = %d\n", jtot, mpi->npy);
    return 1;
  }
  if(ktot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR ktot = %d is not a multiple of npx = %d\n", ktot, mpi->npx);
    return 1;
  }

  imax   = itot / mpi->npx;
  jmax   = jtot / mpi->npy;
  kmax   = ktot;

  iblock = itot / mpi->npy;
  jblock = jtot / mpi->npx;
  kblock = ktot / mpi->npx;

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

  x     = new double[imax+2*igc];
  xh    = new double[imax+2*igc];
  y     = new double[jmax+2*jgc];
  yh    = new double[jmax+2*jgc];
  z     = new double[kmax+2*kgc];
  zh    = new double[kmax+2*kgc];
  dz    = new double[kmax+2*kgc];
  dzh   = new double[kmax+2*kgc];
  dzi   = new double[kmax+2*kgc];
  dzhi  = new double[kmax+2*kgc];
  dzi4  = new double[kmax+2*kgc];
  dzhi4 = new double[kmax+2*kgc];

  allocated = true;

  // initialize the communication functions
  initmpi();

  return 0;
}

int cgrid::create(cinput *inputin)
{
  if(inputin->getProf(&z[kstart], "z", kmax))
    return 1;

  calculate();

  return 0;
}

int cgrid::calculate()
{
  int i,j,k;

  dx = xsize / itot;
  dy = ysize / jtot;

  double xoff = mpi->mpicoordx * xsize / mpi->npx;
  double yoff = mpi->mpicoordy * ysize / mpi->npy;

  // calculate the x and y coordinates
  for(i=0; i<icells; i++)
  {
    x [i] = 0.5*dx + (i-igc)*dx + xoff;
    xh[i] = (i-igc)*dx + xoff;
  }

  for(j=0; j<jcells; j++)
  {
    y [j] = 0.5*dy + (j-jgc)*dy + yoff;
    yh[j] = (j-jgc)*dy + yoff;
  }

  // calculate the height of the ghost cell
  z[kstart-1] = -z[kstart  ];
  z[kstart-2] = -z[kstart+1];
  z[kstart-3] = -z[kstart+2];

  z[kend  ] = 2.*zsize - z[kend-1];
  z[kend+1] = 2.*zsize - z[kend-2];
  z[kend+2] = 2.*zsize - z[kend-3];

  // calculate the half levels according to the numerical scheme
  // zh[0] = -999.;
  // zh[kstart-1] = interp4biasbot(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1]);
  zh[kstart  ] = 0.;
  for(k=kstart+1; k<kend; k++)
    zh[k] = interp4(z[k-2], z[k-1], z[k], z[k+1]);
  zh[kend] = zsize;
  // zh[kend+1] = interp4biastop(z[kend-2], z[kend-1], z[kend], z[kend+1]);
  //
  zh[kstart-1] = -zh[kstart+1];
  zh[kstart-2] = -zh[kstart+2];
  zh[kstart-3] = -zh[kstart+3];

  zh[kend+1] = 2.*zsize - zh[kend-1];
  zh[kend+2] = 2.*zsize - zh[kend-2];

  // compute the height of the grid cells
  for(k=1; k<kcells; k++)
  {
    dzh [k] = z[k] - z[k-1];
    dzhi[k] = 1./dzh[k];
  }
  dzh [kstart-3] = dzh [kstart+3];
  dzhi[kstart-3] = dzhi[kstart+3];

  // compute the height of the grid cells
  for(k=1; k<kcells-1; k++)
  {
    dz [k] = zh[k+1] - zh[k];
    dzi[k] = 1./dz[k];
  }
  dz [kstart-3] = dz [kstart+2];
  dzi[kstart-3] = dzi[kstart+2];
  dz [kend+2] = dz [kend-2];
  dzi[kend+2] = dzi[kend-2];

  /*
  // calculate the inverse gradients for the 4th order scheme
  dzi4 [kstart] = 1./grad4xbiasbot(zh[kstart  ], zh[kstart+1], zh[kstart+2], zh[kstart+3]);
  dzhi4[kstart] = 1./grad4xbiasbot(z [kstart-1], z [kstart  ], z [kstart+1], z [kstart+2]);
  for(k=kstart+1; k<kend-1; k++)
  {
    dzi4 [k] = 1./grad4x(zh[k-1], zh[k  ], zh[k+1], zh[k+2]);
    dzhi4[k] = 1./grad4x(z [k-2], z [k-1], z [k  ], z [k+1]);
  }
  dzi4 [kend-1] = 1./grad4xbiastop(zh[kend-3], zh[kend-2], zh[kend-1], zh[kend]);

  dzhi4[kend-1] = 1./grad4x       (z[kend-3], z[kend-2], z[kend-1], z[kend]);
  dzhi4[kend  ] = 1./grad4xbiastop(z[kend-3], z[kend-2], z[kend-1], z[kend]);
  */
  
  // calculate the inverse gradients for the 4th order scheme
  // dzi4 [0] = -999.;
  // dzhi4[0] = -999.;
  // dzi4 [kstart-1] = 1./grad4xbiasbot(zh[kstart-1], zh[kstart  ], zh[kstart+1], zh[kstart+2]);
  // dzhi4[kstart-1] = 1./grad4xbiasbot(z [kstart-2], z [kstart-1], z [kstart  ], z [kstart+1]);
  //
  for(k=kstart; k<kend; k++)
  {
    dzi4 [k] = 1./grad4x(zh[k-1], zh[k  ], zh[k+1], zh[k+2]);
    dzhi4[k] = 1./grad4x(z [k-2], z [k-1], z [k  ], z [k+1]);
  }
  dzhi4[kend  ] = 1./grad4x(z [kend-2], z [kend-1], z [kend], z [kend+1]);

  // bc's
  dzi4 [kstart-3] = dzi4 [kstart+2];
  dzhi4[kstart-3] = dzhi4[kstart+3];
  dzi4 [kstart-2] = dzi4 [kstart+1];
  dzhi4[kstart-2] = dzhi4[kstart+2];
  dzi4 [kstart-1] = dzi4 [kstart  ];
  dzhi4[kstart-1] = dzhi4[kstart+1];

  dzi4 [kend  ] = dzi4 [kend-1];
  dzi4 [kend+1] = dzi4 [kend-2];
  dzhi4[kend+1] = dzhi4[kend-1];
  dzi4 [kend+2] = dzi4 [kend-3];
  dzhi4[kend+2] = dzhi4[kend-2];

  return 0;
}

inline double cgrid::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

inline double cgrid::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b)); 
}
