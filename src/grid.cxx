#include <cstdio>
#include <cmath>
#include "grid.h"
#include "input.h"

// build the grid
cgrid::cgrid()
{
  std::printf("Creating instance of object grid\n");
  allocated = false;
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
  }

  std::printf("Destroying instance of object grid\n");
}

int cgrid::readinifile(cinput *inputin)
{
  /*// setup Taylor-Green vortex
  xsize = 1.;
  ysize = 1.;
  zsize = 0.5;

  itot  = 64;
  jtot  = 8;
  ktot  = 32;
  // end setup Taylor-Green vortex*/

  int n = 0;

  n += inputin->getItem(&xsize, "grid", "xsize");
  n += inputin->getItem(&ysize, "grid", "ysize");
  n += inputin->getItem(&zsize, "grid", "zsize");

  n += inputin->getItem(&itot, "grid", "itot");
  n += inputin->getItem(&jtot, "grid", "jtot");
  n += inputin->getItem(&ktot, "grid", "ktot");

  if(n > 0)
    return 1;
  
  igc = 1;
  jgc = 1;
  kgc = 1;

  return 0;
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

  x    = new double[imax+2*igc];
  xh   = new double[imax+2*igc];
  y    = new double[jmax+2*jgc];
  yh   = new double[jmax+2*jgc];
  z    = new double[kmax+2*kgc];
  zh   = new double[kmax+2*kgc];
  dz   = new double[kmax+2*kgc];
  dzh  = new double[kmax+2*kgc];
  dzi  = new double[kmax+2*kgc];
  dzhi = new double[kmax+2*kgc];

  allocated = true;

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
  }
  // end Moser180 setup 
  
  /*// uniform height setup
  for(k=kstart; k<kend; k++)
    z[k] = zsize / (2*kmax) + zsize / kmax * (k-kstart);
  // end uniform height setup*/

  calculate();

  return 0;
}

int cgrid::calculate()
{
  int i,j,k;

  dx = xsize / itot;
  dy = ysize / jtot;

  // calculate the x and y coordinates
  for(i=0; i<icells; i++)
  {
    x [i] = 0.5*dx + (i-igc)*dx;
    xh[i] = (i-igc)*dx;
  }

  for(j=0; j<jcells; j++)
  {
    y [j] = 0.5*dy + (j-jgc)*dy;
    yh[j] = (j-jgc)*dy;
  }

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

  // set the non-initialized values
  zh  [0] = -zh[2];
  dzh [0] = -999.;
  dzhi[0] = -999.;

  // compute the heigth of the grid cells
  for(k=kstart; k<kend; k++)
  {
    dz [k] = 0.5*(z[k]-z[k-1]) + 0.5*(z[k+1]-z[k]);
    dzi[k] = 1./dz[k];
  }

  // compute the height of the ghost cells
  for(k=0; k<kgc; k++)
  {
    dz[kstart-k-1]  = dz[kstart+k];
    dz[kend+k]      = dz[kend-k-1];
    dzi[kstart-k-1] = 1./dz[kstart-k-1];
    dzi[kend+k]     = 1./dz[kend+k];
  }

  return 0;
}

int cgrid::save()
{
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" cannot be written\n", filename);
    return 1;
  }
  else
    std::printf("Saving \"%s\"\n", filename);

  fwrite(&x [istart], sizeof(double), itot, pFile);
  fwrite(&xh[istart], sizeof(double), itot, pFile);
  fwrite(&y [jstart], sizeof(double), jtot, pFile);
  fwrite(&yh[jstart], sizeof(double), jtot, pFile);
  fwrite(&z [kstart], sizeof(double), ktot, pFile);
  fwrite(&zh[kstart], sizeof(double), ktot, pFile);
  fclose(pFile);

  return 0;
}

int cgrid::load()
{
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }
  else
    std::printf("Loading \"%s\"\n", filename);

  fread(&x [istart], sizeof(double), itot, pFile);
  fread(&xh[istart], sizeof(double), itot, pFile);
  fread(&y [jstart], sizeof(double), jtot, pFile);
  fread(&yh[jstart], sizeof(double), jtot, pFile);
  fread(&z [kstart], sizeof(double), ktot, pFile);
  fread(&zh[kstart], sizeof(double), ktot, pFile);
  fclose(pFile);

  calculate();

  return 0;
}
