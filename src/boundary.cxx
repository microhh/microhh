#include <cstdio>
#include <cmath>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary.h"
#include "defines.h"

cboundary::cboundary(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object boundary\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cboundary::~cboundary()
{
  // std::printf("Destroying instance of object boundary\n");
}

int cboundary::readinifile(cinput *inputin)
{
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&iboundary, "physics", "iboundary");

  n += inputin->getItem(&bcbotmom, "fields", "bcbotmom");
  n += inputin->getItem(&bctopmom, "fields", "bctopmom");

  n += inputin->getItem(&bcbotscal, "fields", "bcbotscal");
  n += inputin->getItem(&bctopscal, "fields", "bctopscal");

  n += inputin->getItem(&sbot, "fields", "sbot");
  n += inputin->getItem(&stop, "fields", "stop");

  // optional parameters
  n += inputin->getItem(&iboundarytype, "boundary", "iboundarytype", 0);

  // patch type
  n += inputin->getItem(&patch_dim,  "boundary", "patch_dim" , 2 );
  n += inputin->getItem(&patch_xh,   "boundary", "patch_xh"  , 1.);
  n += inputin->getItem(&patch_xr,   "boundary", "patch_xr"  , 1.);
  n += inputin->getItem(&patch_xi,   "boundary", "patch_xi"  , 0.);
  n += inputin->getItem(&patch_facr, "boundary", "patch_facr", 1.);
  n += inputin->getItem(&patch_facl, "boundary", "patch_facl", 0.);

    // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cboundary::setvalues()
{
  setbc((*fields->u).databot, (*fields->u).datagradbot, bcbotmom , 0.  );
  setbc((*fields->v).databot, (*fields->v).datagradbot, bcbotmom , 0.  );

  setbc((*fields->u).datatop, (*fields->u).datagradtop, bctopmom , 0.  );
  setbc((*fields->v).datatop, (*fields->v).datagradtop, bctopmom , 0.  );

  if(iboundarytype == 0)
  {
    setbc((*fields->s).databot, (*fields->s).datagradbot, bcbotscal, sbot);
    setbc((*fields->s).datatop, (*fields->s).datagradtop, bctopscal, stop);
  }
  if(iboundarytype == 1)
  {
    setbc_patch((*fields->s).datagradbot, patch_facl, patch_facr, sbot);
    setbc      ((*fields->s).datatop, (*fields->s).datagradtop, bctopscal, stop);
  }

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

  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij = i + j*jj;
      xmod = fmod(grid->x[i], patch_xh);
      ymod = fmod(grid->y[i], patch_xh);

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

      a[ij] = 0.5*(avall+avalr) + (avalr-avall)*errvalx*errvaly;
    }

  return 0;
}


int cboundary::setbc(double * restrict a, double * restrict agrad, int sw, double aval)
{
  int ij,jj;
  jj = grid->icells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij = i + j*jj;
        a[ij] = aval;
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij = i + j*jj;
        agrad[ij] = aval;
      }
  }

  return 0;
}

int cboundary::exec()
{
  if(iboundary == 2)
  {
    // bottom boundary conditions
    setgcbot_2nd((*fields->u).data, grid->dzh, bcbotmom, 0.);
    setgcbot_2nd((*fields->v).data, grid->dzh, bcbotmom, 0.);

    // top boundary conditions
    setgctop_2nd((*fields->u).data, grid->dzh, bctopmom, 0.);
    setgctop_2nd((*fields->v).data, grid->dzh, bctopmom, 0.);
    setgctop_2nd((*fields->s).data, grid->dzh, bctopscal, stop);
  }
  else if(iboundary == 4)
  {
    // bottom boundary conditions
    setgcbot_4th ((*fields->u).data, grid->z, bcbotmom, 0.);
    setgcbot_4th ((*fields->v).data, grid->z, bcbotmom, 0.);
    //setgcbot_4th ((*fields->s).data, grid->z, bcbotscal, sbot);
    setgcbot_4th ((*fields->s).data, grid->z, bcbotscal, (*fields->s).datagradbot);
    setgcbotw_4th((*fields->w).data);

    // top boundary conditions
    setgctop_4th((*fields->u).data, grid->z, bctopmom, 0.);
    setgctop_4th((*fields->v).data, grid->z, bctopmom, 0.);
    //setgctop_4th((*fields->s).data, grid->z, bctopscal, stop);
    setgctop_4th ((*fields->s).data, grid->z, bctopscal, (*fields->s).datagradtop);
    setgctopw_4th((*fields->w).data);
  }
 
  // cyclic boundary conditions
  grid->boundary_cyclic((*fields->u).data);
  grid->boundary_cyclic((*fields->v).data);
  grid->boundary_cyclic((*fields->w).data);
  grid->boundary_cyclic((*fields->s).data);

  return 0;
}

// BOUNDARY CONDITIONS THAT HAVE ONE VALUE FOR THE ENTIRE DOMAIN
int cboundary::setgcbot_2nd(double * restrict a, double * restrict dzh, int sw, double abot)
{ 
  int ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = 2.*abot - a[ijk];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
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
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        // add the bcvalues later
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = 2.*atop - a[ijk];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        // add the bcvalues later
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
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = (8./3.)*abot - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
        a[ijk-kk2] = 8.*abot - 9.*a[ijk] + 2.*a[ijk+kk1];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
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
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (8./3.)*atop - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
        a[ijk+kk2] = 8.*atop - 9.*a[ijk] + 2.*a[ijk-kk1];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*atop + a[ijk    ];
        a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*atop + a[ijk-kk1];
      }
  }

  return 0;
}

// BOUNDARY CONDITIONS THAT CONTAIN A 2D PATTERN
int cboundary::setgcbot_2nd(double * restrict a, double * restrict dzh, int sw, double * restrict abot)
{ 
  int ij,ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = 2.*abot[ij] - a[ijk];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = -abot[ij]*dzh[kstart] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgctop_2nd(double * restrict a, double * restrict dzh, int sw, double * restrict atop)
{ 
  int ij,ijk,jj,kk,kend;

  kend = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = 2.*atop[ij] - a[ijk];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = atop[ij]*dzh[kend] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgcbot_4th(double * restrict a, double * restrict z, int sw, double * restrict abot)
{ 
  int ij,ijk,jj,kk1,kk2,kstart;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = (8./3.)*abot[ij] - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
        a[ijk-kk2] = 8.*abot[ij] - 9.*a[ijk] + 2.*a[ijk+kk1];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = -(1./24.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*abot[ij] + a[ijk    ];
        a[ijk-kk2] = -(1./ 8.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*abot[ij] + a[ijk+kk1];
      }
  }

  return 0;
}

int cboundary::setgctop_4th(double * restrict a, double * restrict z, int sw, double * restrict atop)
{ 
  int ij,ijk,jj,kend,kk1,kk2;

  kend = grid->kend;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (8./3.)*atop[ij] - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
        a[ijk+kk2] = 8.*atop[ij] - 9.*a[ijk] + 2.*a[ijk-kk1];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*atop[ij] + a[ijk    ];
        a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*atop[ij] + a[ijk-kk1];
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

  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
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

  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
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
