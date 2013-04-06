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

  // read the options for the patch type
  if(swboundarytype == "1")
  {
    // patch type
    n += inputin->getItem(&patch_dim,  "boundary", "patch_dim" , "", 2 );
    n += inputin->getItem(&patch_xh,   "boundary", "patch_xh"  , "", 1.);
    n += inputin->getItem(&patch_xr,   "boundary", "patch_xr"  , "", 1.);
    n += inputin->getItem(&patch_xi,   "boundary", "patch_xi"  , "", 0.);
    n += inputin->getItem(&patch_facr, "boundary", "patch_facr", "", 1.);
    n += inputin->getItem(&patch_facl, "boundary", "patch_facl", "", 0.);
  }
  else if(swboundarytype == "surface")
  {
    n += inputin->getItem(&z0m, "surface", "z0m", "");
    n += inputin->getItem(&z0h, "surface", "z0h", "");
  }

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

// TODO make this the surface init function
int cboundary::init()
{
  if(swboundarytype == "surface")
  {
    obuk  = new double[grid->icells*grid->jcells];
    ustar = new double[grid->icells*grid->jcells];
  }

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
    else if(swboundarytype == "1")
    {
      setbc_patch(it->second->datagradbot, patch_facl, patch_facr, sbc[it->first]->bot);
      setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc);
    }
    // TODO move to surface
    else if(swboundarytype == "surface")
    {
      setbc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc);
      setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc);

      // TODO temporary assignment of ustar, to yield z/L of -z/Ltemp
      int ij,jj;
      jj = grid->icells;

      double Ltemp = -10.;

      for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for(int i=grid->istart; i<grid->iend; ++i)
        {
          ij = i + j*jj;
          ustar[ij] = std::pow(-grid->z[grid->kstart]*kappa*9.81/300.*fields->sp["s"]->datafluxbot[ij]/Ltemp, 1./3.);
        }
    }
  }

  return 0;
}

int cboundary::exec()
{
  if(swboundarytype == "surface")
    surface();

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

// TODO move all underlying functions to their own file
// surface model
int cboundary::surface()
{
  // TODO replace this later to the right location
  // start with retrieving the stability information
  stability(obuk, ustar, fields->sp["s"]->datafluxbot,
            fields->u->data, fields->v->data, fields->sp["s"]->data,
            fields->u->databot, fields->v->databot, fields->sp["s"]->databot,
            fields->s["tmp1"]->data, fields->s["tmp1"]->data,
            grid->z);

  // calculate the surface value, gradient and flux depending on the chosen boundary condition
  surfvalues(ustar, obuk, fields->sp["s"]->data,
             fields->sp["s"]->databot, fields->sp["s"]->datagradbot, fields->sp["s"]->datafluxbot,
             grid->z[grid->kstart]);

  return 0;
}

int cboundary::stability(double * restrict obuk, double * restrict ustar, double * restrict bfluxbot,
                         double * restrict u   , double * restrict v    , double * restrict b, 
                         double * restrict ubot, double * restrict vbot , double * restrict bbot,
                         double * restrict utot, double * restrict ubottot, double * restrict z)
{
  int ij,ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  // calculate total wind
  // first, interpolate the wind to the scalar location
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      ubottot[ij] = std::pow(std::pow(ubot[ij ], 2.) + std::pow(vbot[ij ], 2.), 0.5);
      utot   [ij] = std::pow(std::pow(u   [ijk], 2.) + std::pow(v   [ijk], 2.), 0.5);
    }

  // TODO replace by value from boundary
  double gravitybeta = 9.81/300.;

  // calculate Obukhov length
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      obuk[ij] = -std::pow(ustar[ij], 3.) / (gravitybeta*kappa*bfluxbot[ij]);
    }

  return 0;
}

int cboundary::surfvalues(double * restrict ustar, double * restrict obuk, double * restrict var,
                          double * restrict varbot, double * restrict vargradbot, double * restrict varfluxbot, 
                          double zsl)
{
  int ij,ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  // first calculate the surface value
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      varbot[ij]     = varfluxbot[ij] / (ustar[ij]*fh(zsl, z0h, obuk[ij])) + var[ijk];
      vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
      //if(i==3 && j==3) std::printf("CvH (ustar, bot, fluxbot, gradbot): %d, %d, %E, %E, %E, %E\n", i, j, ustar[ij], varbot[ij], varfluxbot[ij], vargradbot[ij]);
    }

  return 0;
}

/*
double cboundary::calcobuk(double ustar, double wtheta, double du, double db, double zsl, 
                           bool doustar, bool dofstar)
{
  double L, L0;
  double Lstart, Lend;
  double fx, fxdif;

  if(doustar == false && dofstar == false)
    L = -std::pow(ustar, 3.) / (gravitybeta*kappa*zsl*wtheta)

  if(Rib > 0.)
  {
    L  = 1.;
    L0 = 2.;
  }
  else
  {
    L  = -1.;
    L0 = -2.;
  }

  while(std::abs(L - L0) > 0.001)
  {
    L0     = L;
    fx     = Rib - zsl/L * (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L)) / std::pow(std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L), 2.);
    Lstart = L - 0.001*L;
    Lend   = L + 0.001*L;
    fxdif  = ( (- zsl/Lstart * (std::log(zsl/z0h) - psih(zsl/Lstart) + psih(z0h/Lstart)) / pow(std::log(zsl/z0m) - psim(zsl/Lstart) + psim(z0m/Lstart), 2.)) - (-zsl/Lend * (std::log(zsl/z0h) - psih(zsl/Lend) + psih(z0h/Lend)) / pow(std::log(zsl/z0m) - psim(zsl/Lend) + psim(z0m/Lend), 2.)) ) / (Lstart-Lend);
    L      = L - fx/fxdif;
  }

  return L;
}
*/

inline double cboundary::fm(double zsl, double z0m, double L)
{
  double fm;
  fm = kappa / (std::log(zsl/z0m) - psim(zsl/L) + psih(z0m/L));
  return fm;
}

inline double cboundary::fh(double zsl, double z0h, double L)
{
  double fh;
  fh = kappa / (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L));
  return fh;
}

inline double cboundary::psim(double zeta)
{
  double psim;
  double x;
  if(zeta <= 0.)
  {
    // Businger-Dyer functions
    // x     = (1. - 16. * zeta) ** (0.25)
    // psim  = 3.14159265 / 2. - 2. * arctan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
    // Wilson functions
    x    = std::pow(1. + std::pow(3.6 * std::abs(zeta),2./3.), -0.5);
    psim = 3.*std::log( (1. + 1./x) / 2.);
  }
  else
  {
    psim = -2./3.*(zeta - 5./0.35) * exp(-0.35 * zeta) - zeta - (10./3.) / 0.35;
  }
  return psim;
}

inline double cboundary::psih(double zeta)
{
  double psih;
  double x;
  if(zeta <= 0.)
  {
    // Businger-Dyer functions
    // x     = (1. - 16. * zeta) ** (0.25)
    // psih  = 2. * log( (1. + x ** 2.) / 2. )
    // Wilson functions
    x    = std::pow(1. + std::pow(7.9*std::abs(zeta), (2./3.)), -0.5);
    psih = 3. * std::log( (1. + 1. / x) / 2.);
  }
  else
  {
    psih  = -2./3. * (zeta-5./0.35) * exp(-0.35*zeta) - std::pow(1. + (2./3.) * zeta, 1.5) - (10./3.) / 0.35 + 1.;
  }
  return psih;
}

inline double cboundary::phim(double zeta)
{
  double phim;
  if(zeta <= 0.)
  {
    // Businger-Dyer functions
    //x     = (1. - 16. * zeta) ** (0.25)
    //psim  = 3.14159265 / 2. - 2. * arctan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
    // Wilson functions
    phim = std::pow(1. + 3.6*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }
  else
    phim = 1. + 5.*zeta;

  return phim;
}

inline double cboundary::phih(double zeta)
{
  double phih;
  if(zeta <= 0.)
  {
    // Businger-Dyer functions
    // x     = (1. - 16. * zeta) ** (0.25)
    // psih  = 2. * log( (1. + x ** 2.) / 2. )
    // Wilson functions
    phih = std::pow(1. + 7.9*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }
  else
    phih = 1. + 5.*zeta;

  return phih;
}
