#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "diff.h"
#include "defines.h"

cdiff::cdiff(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object diff\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  diff_g2  = new cdiff_g2 (grid, fields, mpi);
  diff_g42 = new cdiff_g42(grid, fields, mpi);
  diff_g4  = new cdiff_g4 (grid, fields, mpi);
}

cdiff::~cdiff()
{
  delete diff_g2;
  delete diff_g42;
  delete diff_g4;

  std::printf("Destroying instance of object diff\n");
}

int cdiff::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&idiff, "physics", "idiff");

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cdiff::init()
{
  // get the maximum time step for diffusion
  double viscmax = std::max(fields->visc, fields->viscs);

  dnmul = 0;
  for(int k=grid->kstart; k<grid->kend; k++)
    dnmul = std::max(dnmul, std::abs(viscmax * (1./(grid->dx*grid->dx) + 1./(grid->dy*grid->dy) + 1./(grid->dz[k]*grid->dz[k]))));

  return 0;
}

double cdiff::getdn(double dt)
{
  double dn;

  // in case of no diffusion set dn to a small number to avoid zero divisions
  if(idiff == 0)
    dn = dsmall;
  else
    dn = dnmul*dt;

  return dn;
}

int cdiff::exec()
{
  if(idiff == 0)
    return 0;

  // diffuse the flow
  if(idiff == 2)
  {
    diff_g2->diffc((*fields->ut).data, (*fields->u).data, grid->dzi, grid->dzhi, fields->visc);
    diff_g2->diffc((*fields->vt).data, (*fields->v).data, grid->dzi, grid->dzhi, fields->visc);
    diff_g2->diffw((*fields->wt).data, (*fields->w).data, grid->dzi, grid->dzhi, fields->visc);

    diff_g2->diffc((*fields->st).data, (*fields->s).data, grid->dzi, grid->dzhi, fields->viscs);
  }
  else if(idiff == 42)
  {
    diff_g42->diffc((*fields->ut).data, (*fields->u).data, grid->dzi, grid->dzhi, fields->visc);
    diff_g42->diffc((*fields->vt).data, (*fields->v).data, grid->dzi, grid->dzhi, fields->visc);
    diff_g42->diffw((*fields->wt).data, (*fields->w).data, grid->dzi, grid->dzhi, fields->visc);

    diff_g42->diffc((*fields->st).data, (*fields->s).data, grid->dzi, grid->dzhi, fields->viscs);
  }
  else if(idiff == 4)
  {
    diff_g4->diffc((*fields->ut).data, (*fields->u).data, grid->z, grid->zh, fields->visc);
    diff_g4->diffc((*fields->vt).data, (*fields->v).data, grid->z, grid->zh, fields->visc);
    diff_g4->diffw((*fields->wt).data, (*fields->w).data, grid->dzi, grid->dzhi, fields->visc);

    diff_g4->diffc((*fields->st).data, (*fields->s).data, grid->z, grid->zh, fields->viscs);
  }

  return 0;
}

