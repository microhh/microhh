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
  // std::printf("Creating instance of object diff\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  diff_g2  = new cdiff_g2 (grid, fields, mpi);
  diff_g42 = new cdiff_g42(grid, fields, mpi);
  diff_g4  = new cdiff_g4 (grid, fields, mpi);

  diff_les_g2  = new cdiff_les_g2(grid, fields, mpi);
}

cdiff::~cdiff()
{
  delete diff_g2;
  delete diff_g42;
  delete diff_g4;

  delete diff_les_g2;
  // std::printf("Destroying instance of object diff\n");
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

int cdiff::setvalues()
{
  // get the maximum time step for diffusion
  double viscmax = fields->visc;
  for(fieldmap::iterator it = fields->sp.begin(); it!=fields->sp.end(); it++)
    viscmax = std::max(it->second->visc, viscmax);

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
  if(idiff == 22)
  {
    // calculate eddy viscosity
    diff_les_g2->evisc((*fields->evisc).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, (*fields->s["s"]).data, grid->z, grid->dz, grid->dzi, grid->dzhi, fields->tPr);
    dn = diff_les_g2->getdn((*fields->evisc).data, grid->dzi, fields->tPr);
  }
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

    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      diff_g2->diffc((*it->second).data, (*fields->s[it->first]).data, grid->dzi, grid->dzhi, (*it->second).visc);
  }
  else if(idiff == 42)
  {
    diff_g42->diffc((*fields->ut).data, (*fields->u).data, grid->dzi, grid->dzhi, fields->visc);
    diff_g42->diffc((*fields->vt).data, (*fields->v).data, grid->dzi, grid->dzhi, fields->visc);
    diff_g42->diffw((*fields->wt).data, (*fields->w).data, grid->dzi, grid->dzhi, fields->visc);


    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      diff_g42->diffc((*it->second).data, (*fields->s[it->first]).data, grid->dzi, grid->dzhi, (*it->second).visc);
  }
  else if(idiff == 4)
  {
    diff_g4->diffc((*fields->ut).data, (*fields->u).data, grid->dzi4, grid->dzhi4, fields->visc);
    diff_g4->diffc((*fields->vt).data, (*fields->v).data, grid->dzi4, grid->dzhi4, fields->visc);
    diff_g4->diffw((*fields->wt).data, (*fields->w).data, grid->dzi4, grid->dzhi4, fields->visc);


    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      diff_g4->diffc((*it->second).data, (*fields->s[it->first]).data, grid->dzi4, grid->dzhi4, (*it->second).visc);
  }

  else if(idiff == 22)
  {
    diff_les_g2->diffu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, grid->dzhi, (*fields->evisc).data, (*fields->u).datafluxbot, (*fields->u).datafluxtop);
    diff_les_g2->diffv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, grid->dzhi, (*fields->evisc).data, (*fields->v).datafluxbot, (*fields->v).datafluxtop);
    diff_les_g2->diffw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, grid->dzhi, (*fields->evisc).data);

    diff_les_g2->diffc((*fields->st["s"]).data, (*fields->s["s"]).data, grid->dzi, grid->dzhi, (*fields->evisc).data, (*fields->s["s"]).datafluxbot, (*fields->s["s"]).datafluxtop, fields->tPr);
  }

  return 0;
}

