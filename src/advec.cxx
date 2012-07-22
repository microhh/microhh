#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "defines.h"

cadvec::cadvec(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object advec\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  advec_g2   = new cadvec_g2  (grid, fields, mpi);
  advec_g2i4 = new cadvec_g2i4(grid, fields, mpi);
}

cadvec::~cadvec()
{
  delete advec_g2;
  delete advec_g2i4;

  std::printf("Destroying instance of object advec\n");
}

int cadvec::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&iadvec, "physics", "iadvec");

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

double cadvec::getcfl(double dt)
{
  double cfl;

  // in case of no advection set cfl to a small number to avoid zero divisions
  if(iadvec == 0)
    cfl = dsmall;
  
  if(iadvec == 2)
    cfl = advec_g2->calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);
  else if(iadvec == 4)
    cfl = advec_g2i4->calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);

  return cfl;
}

int cadvec::exec()
{
  if(iadvec == 0)
    return 0;

  if(iadvec == 2)
  {
    advec_g2->advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2->advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2->advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi);
    advec_g2->advecs((*fields->st).data, (*fields->s).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
  }
  else if(iadvec == 4)
  {
    advec_g2i4->advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2i4->advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2i4->advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi);
    advec_g2i4->advecs((*fields->st).data, (*fields->s).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
  }

  return 0;
}

