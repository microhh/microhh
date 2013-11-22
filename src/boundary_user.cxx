#include <cstdio>
#include <cmath>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_user.h"
#include "defines.h"

#define NO_VELOCITY 0.
#define NO_OFFSET 0.

cboundary_user::cboundary_user(cgrid *gridin, cfields *fieldsin, cmpi *mpiin) : cboundary(gridin, fieldsin, mpiin)
{
}

int cboundary_user::readinifile(cinput *inputin)
{
  int n = 0;

  // obligatory parameters
  // n += inputin->getItem(&swboundary, "boundary", "swboundary", "", grid->swspatialorder);
  swspatialorder = grid->swspatialorder;

  n += inputin->getItem(&mbcbot, "boundary", "mbcbot", "");
  n += inputin->getItem(&mbctop, "boundary", "mbctop", "");

  // read the boundaries per field
  for(fieldmap::iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    sbc[it->first] = new field3dbc;
    n += inputin->getItem(&sbc[it->first]->bcbot, "boundary", "sbcbot", it->first);
    n += inputin->getItem(&sbc[it->first]->bctop, "boundary", "sbctop", it->first);
    n += inputin->getItem(&sbc[it->first]->bot  , "boundary", "sbot"  , it->first);
    n += inputin->getItem(&sbc[it->first]->top  , "boundary", "stop"  , it->first);
  }

  // patch type
  n += inputin->getItem(&patch_dim,  "boundary", "patch_dim" , "", 2 );
  n += inputin->getItem(&patch_xh,   "boundary", "patch_xh"  , "", 1.);
  n += inputin->getItem(&patch_xr,   "boundary", "patch_xr"  , "", 1.);
  n += inputin->getItem(&patch_xi,   "boundary", "patch_xi"  , "", 0.);
  n += inputin->getItem(&patch_facr, "boundary", "patch_facr", "", 1.);
  n += inputin->getItem(&patch_facl, "boundary", "patch_facl", "", 0.);
  
  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cboundary_user::setvalues()
{
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->u);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->v);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->u);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->v);

  for(fieldmap::iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setbc_patch(it->second->datagradbot, patch_facl, patch_facr, sbc[it->first]->bot);
    setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, NO_OFFSET);
  }

  return 0;
}

int cboundary_user::setbc_patch(double * restrict a, double facl, double facr, double aval)
{
  double avall, avalr;
  double xmod, ymod;
  double errvalx, errvaly;

  int ij,jj;
  jj = grid->icells;

  avall = facl*aval;
  avalr = facr*aval;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij = i + j*jj;
      xmod = fmod(grid->x[i], patch_xh);
      ymod = fmod(grid->y[j], patch_xh);

      errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

      if(patch_dim == 2)
        errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);
      else
        errvaly = 1.;

      a[ij] = avall + (avalr-avall)*errvalx*errvaly;
    }

  return 0;
}
