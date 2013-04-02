#include <cstdio>
#include <cmath>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary.h"
#include "defines.h"

cboundary::cboundary(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cboundary::~cboundary()
{
  for(bcmap::iterator it=sbc.begin(); it!=sbc.end(); ++it)
    delete it->second;

  // empty the map
  sbc.clear();
}

int cboundary::readinifile(cinput *inputin)
{
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&swboundary, "boundary", "swboundary", "", grid->swspatialorder);

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

  // optional parameters
  n += inputin->getItem(&swboundarytype, "boundary", "swboundarytype", "", "0");

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

int cboundary::setvalues()
{
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, 0., fields->visc);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, 0., fields->visc);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, 0., fields->visc);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, 0., fields->visc);

  for(fieldmap::iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    if(swboundarytype == "0")
    {
      setbc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc);
      setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc);
    }
    if(swboundarytype == "1")
    {
      setbc_patch(it->second->datagradbot, patch_facl, patch_facr, sbc[it->first]->bot);
      setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc);
    }
  }
  return 0;
}

int cboundary::exec()
{
  if(swboundary == "2")
  {
    setgcbot_2nd(fields->u->data, grid->dzh, mbcbot, fields->u->databot, fields->u->datagradbot);
    setgctop_2nd(fields->u->data, grid->dzh, mbctop, fields->u->datatop, fields->u->datagradtop);

    setgcbot_2nd(fields->v->data, grid->dzh, mbcbot, fields->v->databot, fields->v->datagradbot);
    setgctop_2nd(fields->v->data, grid->dzh, mbctop, fields->v->datatop, fields->v->datagradtop);
    
    for(fieldmap::iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      setgcbot_2nd(it->second->data, grid->dzh, sbc[it->first]->bcbot, it->second->databot, it->second->datagradbot);
      setgctop_2nd(it->second->data, grid->dzh, sbc[it->first]->bctop, it->second->datatop, it->second->datagradtop);
    }
  }
  else if(swboundary == "4")
  {
    setgcbot_4th (fields->u->data, grid->z, mbcbot, fields->u->databot, fields->u->datagradbot);
    setgctop_4th (fields->u->data, grid->z, mbctop, fields->u->datatop, fields->u->datagradtop);

    setgcbot_4th (fields->v->data, grid->z, mbcbot, fields->v->databot, fields->v->datagradbot);
    setgctop_4th (fields->v->data, grid->z, mbctop, fields->v->datatop, fields->v->datagradtop);

    setgcbotw_4th(fields->w->data);
    setgctopw_4th(fields->w->data);
    
    for(fieldmap::iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      setgcbot_4th(it->second->data, grid->z, sbc[it->first]->bcbot, it->second->databot, it->second->datagradbot);
      setgctop_4th(it->second->data, grid->z, sbc[it->first]->bctop, it->second->datatop, it->second->datagradtop);
    }
  }

  // cyclic boundary conditions
  grid->boundary_cyclic(fields->u->data);
  grid->boundary_cyclic(fields->v->data);
  grid->boundary_cyclic(fields->w->data);
  
  for(fieldmap::iterator it = fields->sp.begin(); it!=fields->sp.end(); ++it)
    grid->boundary_cyclic(it->second->data);

  return 0;
}

int cboundary::setbc_patch(double * restrict a, double facl, double facr, double aval)
{
  // dimensions patches
  double xrmid   = 0.5*patch_xh;
  double xrstart = 0.5*(patch_xh - patch_xr);
  double xrend   = 0.5*(patch_xh + patch_xr);

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

      if(xmod < xrmid)
        errvalx =  0.5*erf(0.5*(xmod-xrstart) / patch_xi);
      else
        errvalx = -0.5*erf(0.5*(xmod-xrend) / patch_xi);

      if(patch_dim == 2)
      {
        if(ymod < xrmid)
          errvaly =  0.5*erf(0.5*(ymod-xrstart) / patch_xi);
        else
          errvaly = -0.5*erf(0.5*(ymod-xrend) / patch_xi);
      }
      else
        errvaly = 1.;

      // normalize the values between 0 and 1
      errvalx = errvalx + 0.5;
      errvaly = errvaly + 0.5;

      a[ij] = avall + (avalr-avall)*errvalx*errvaly;
    }

  return 0;
}

int cboundary::setbc(double * restrict a, double * restrict agrad, double * restrict aflux, int sw, double aval, double visc)
{
  int ij,jj;
  jj = grid->icells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        a[ij] = aval;
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        agrad[ij] = aval;
        aflux[ij] = -aval*visc;
      }
  }
  else if(sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        aflux[ij] = aval;
        agrad[ij] = -aval/visc;
      }
  }

  return 0;
}

/*
// BOUNDARY CONDITIONS THAT HAVE ONE VALUE FOR THE ENTIRE DOMAIN
int cboundary::setgcbot_2nd(double * restrict a, double * restrict dzh, int sw, double abot)
{ 
  int ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = 2.*abot - a[ijk];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = -abot*dzh[kstart] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgctop_2nd(double * restrict a, double * restrict dzh, int sw, double atop)
{ 
  int ijk,jj,kk,kend;

  kend = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = 2.*atop - a[ijk];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = atop*dzh[kend] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgcbot_4th(double * restrict a, double * restrict z, int sw, double abot)
{ 
  int ijk,jj,kk1,kk2,kstart;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = (8./3.)*abot - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
        a[ijk-kk2] = 8.*abot - 9.*a[ijk] + 2.*a[ijk+kk1];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = -(1./24.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*abot + a[ijk    ];
        a[ijk-kk2] = -(1./ 8.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*abot + a[ijk+kk1];
      }
  }

  return 0;
}

int cboundary::setgctop_4th(double * restrict a, double * restrict z, int sw, double atop)
{ 
  int ijk,jj,kend,kk1,kk2;

  kend = grid->kend;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (8./3.)*atop - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
        a[ijk+kk2] = 8.*atop - 9.*a[ijk] + 2.*a[ijk-kk1];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*atop + a[ijk    ];
        a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*atop + a[ijk-kk1];
      }
  }

  return 0;
}
*/

// BOUNDARY CONDITIONS THAT CONTAIN A 2D PATTERN
int cboundary::setgcbot_2nd(double * restrict a, double * restrict dzh, int sw, double * restrict abot, double * restrict agradbot)
{ 
  int ij,ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = 2.*abot[ij] - a[ijk];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgctop_2nd(double * restrict a, double * restrict dzh, int sw, double * restrict atop, double * restrict agradtop)
{ 
  int ij,ijk,jj,kk,kend;

  kend = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = 2.*atop[ij] - a[ijk];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgcbot_4th(double * restrict a, double * restrict z, int sw, double * restrict abot, double * restrict agradbot)
{ 
  int ij,ijk,jj,kk1,kk2,kstart;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = (8./3.)*abot[ij] - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
        a[ijk-kk2] = 8.*abot[ij] - 9.*a[ijk] + 2.*a[ijk+kk1];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = -(1./24.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
        a[ijk-kk2] = -(1./ 8.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
      }
  }

  return 0;
}

int cboundary::setgctop_4th(double * restrict a, double * restrict z, int sw, double * restrict atop, double * restrict agradtop)
{ 
  int ij,ijk,jj,kend,kk1,kk2;

  kend = grid->kend;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (8./3.)*atop[ij] - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
        a[ijk+kk2] = 8.*atop[ij] - 9.*a[ijk] + 2.*a[ijk-kk1];
      }
  }
  else if(sw == 1 || sw == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
        a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
      }
  }

  return 0;
}

// BOUNDARY CONDITIONS FOR THE VERTICAL VELOCITY (NO PENETRATION)
int cboundary::setgcbotw_4th(double * restrict w)
{ 
  int ijk,jj,kk1,kk2,kstart;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kstart = grid->kstart;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ijk = i + j*jj + kstart*kk1;
      w[ijk-kk1] = -w[ijk+kk1];
      w[ijk-kk2] = -w[ijk+kk2];
    }
 
  return 0;
}

int cboundary::setgctopw_4th(double * restrict w)
{ 
  int ijk,jj,kk1,kk2,kend;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kend = grid->kend;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ijk = i + j*jj + kend*kk1;
      w[ijk+kk1] = -w[ijk-kk1];
      w[ijk+kk2] = -w[ijk-kk2];
    }
 
  return 0;
}

inline double cboundary::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b)); 
}
