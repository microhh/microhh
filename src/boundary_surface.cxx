#include <cstdio>
#include <cmath>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_surface.h"
#include "defines.h"

// a sign function
inline double sign(double n) { return n > 0 ? 1 : (n < 0 ? -1 : 0);}

cboundary_surface::cboundary_surface(cgrid *gridin, cfields *fieldsin, cmpi *mpiin) : cboundary(gridin, fieldsin, mpiin)
{
  allocated = false;
}

cboundary_surface::~cboundary_surface()
{
  if(allocated)
  {
    delete[] ustar;
    delete[] obuk;
  }
}

int cboundary_surface::readinifile(cinput *inputin)
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

  n += inputin->getItem(&z0m, "boundary", "z0m", "");
  n += inputin->getItem(&z0h, "boundary", "z0h", "");

  // copy all the boundary options and set the model ones to flux type
  surfmbcbot = mbcbot;
  mbcbot = 2;

  // crash in case fixed gradient is prescribed
  if(surfmbcbot == 1)
  {
    if(mpi->mpiid == 0) std::printf("ERROR fixed gradient bc is not supported in surface model\n");
    ++n;
  }
  // read the ustar value only if fixed fluxes are prescribed
  else if(surfmbcbot == 2)
    n += inputin->getItem(&ustarin, "boundary", "ustar", "");

  // process the scalars
  for(bcmap::iterator it=sbc.begin(); it!=sbc.end(); ++it)
  {
    surfsbcbot[it->first] = it->second->bcbot;
    it->second->bcbot = 2;

    // crash in case fixed gradient is prescribed
    if(surfsbcbot[it->first] == 1)
    {
      if(mpi->mpiid == 0) std::printf("ERROR fixed gradient bc is not supported in surface model\n");
      ++n;
    }

    // crash in case of fixed momentum flux and dirichlet bc for scalar
    if(surfsbcbot[it->first] == 0 && surfmbcbot == 2)
    {
      if(mpi->mpiid == 0) std::printf("ERROR fixed ustar bc in combination with Dirichlet bc for scalars is not supported\n");
      ++n;
    }
    // crash in case of no-slip bc for momentum and dirichlet bc for scalar
    else if(surfsbcbot[it->first] == 0 && surfmbcbot == 0)
    {
      if(mpi->mpiid == 0) std::printf("ERROR no-slip bc in combination with Dirichlet bc for scalars is not supported\n");
      ++n;
    }
  }

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cboundary_surface::init()
{
  obuk  = new double[grid->icells*grid->jcells];
  ustar = new double[grid->icells*grid->jcells];

  allocated = true;

  int ij,jj;
  jj = grid->icells;

  // initialize the obukhov length on a small number
  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
   for(int i=0; i<grid->icells; ++i)
   {
     ij = i + j*jj;
     obuk[ij] = dsmall;
   }

  return 0;
}

int cboundary_surface::save(int iotime)
{
  char filename[256];

  std::sprintf(filename, "obuk.%07d", iotime);
  if(mpi->mpiid == 0) std::printf("Saving \"%s\"\n", filename);
  if(grid->savexyslice(obuk, fields->s["tmp1"]->data, filename))
    return 1;

  return 0;
}

int cboundary_surface::load(int iotime)
{
  char filename[256];

  std::sprintf(filename, "obuk.%07d", iotime);
  if(mpi->mpiid == 0) std::printf("Loading \"%s\"\n", filename);
  if(grid->loadxyslice(obuk, fields->s["tmp1"]->data, filename))
    return 1;

  grid->boundary_cyclic2d(obuk);

  return 0;
}

int cboundary_surface::setvalues()
{
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, 0., fields->visc);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, 0., fields->visc);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, 0., fields->visc);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, 0., fields->visc);

  for(fieldmap::iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setbc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc);
    setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc);
  }

  int ij,jj;
  jj = grid->icells;

  // in case the momentum has a fixed ustar, set the value to that of the input
  if(surfmbcbot == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
        {
          ij = i + j*jj;
          // limit ustar at 1e-4 to avoid zero divisions
          ustar[ij] = std::max(0.0001, ustarin);
        }
   }

  return 0;
}

// surface model
int cboundary_surface::bcvalues()
{
  // start with retrieving the stability information
  // TODO make this working properly with buoyancy
  stability(ustar, obuk, fields->sp["s"]->datafluxbot,
            fields->u->data, fields->v->data, fields->sp["s"]->data,
            fields->u->databot, fields->v->databot, fields->sp["s"]->databot,
            fields->s["tmp1"]->data, grid->z);

  // calculate the surface value, gradient and flux depending on the chosen boundary condition
  surfm(ustar, obuk,
        fields->u->data, fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot,
        fields->v->data, fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot,
        grid->z[grid->kstart], surfmbcbot);

  surfs(ustar, obuk, fields->sp["s"]->data,
        fields->sp["s"]->databot, fields->sp["s"]->datagradbot, fields->sp["s"]->datafluxbot,
        grid->z[grid->kstart], surfsbcbot["s"]);

  return 0;
}

int cboundary_surface::stability(double * restrict ustar, double * restrict obuk, double * restrict bfluxbot,
                                 double * restrict u    , double * restrict v   , double * restrict b       ,
                                 double * restrict ubot , double * restrict vbot, double * restrict bbot    ,
                                 double * restrict dutot, double * restrict z)
{
  int ij,ijk,ii,jj,kk,kstart;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  // calculate total wind
  double utot, ubottot;
  const double minval = 1.e-1;
  // first, interpolate the wind to the scalar location
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      ubottot = std::pow(  0.5*(std::pow(ubot[ij], 2.) + std::pow(ubot[ij+ii], 2.))
                         + 0.5*(std::pow(vbot[ij], 2.) + std::pow(vbot[ij+jj], 2.)), 0.5);
      utot    = std::pow(  0.5*(std::pow(u[ijk], 2.) + std::pow(u[ijk+ii], 2.))
                         + 0.5*(std::pow(v[ijk], 2.) + std::pow(v[ijk+jj], 2.)), 0.5);
      // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
      // otherwise evisc at k = kstart blows up
      dutot[ij] = std::max(std::abs(utot - ubottot), minval);
    }

  grid->boundary_cyclic2d(dutot);

  // TODO replace by value from buoyancy
  double gravitybeta = 9.81/300.;

  // calculate Obukhov length
  // case 1: fixed buoyancy flux and fixed ustar
  if(surfmbcbot == 2 && surfsbcbot["s"] == 2)
  {
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        obuk[ij] = -std::pow(ustar[ij], 3.) / (gravitybeta*kappa*bfluxbot[ij]);
      }
  }
  // case 2: fixed buoyancy surface value and free ustar
  if(surfmbcbot == 0 && surfsbcbot["s"] == 2)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        obuk [ij] = calcobuk(obuk[ij], dutot[ij], bfluxbot[ij], z[kstart]);
        ustar[ij] = dutot[ij] * fm(z[kstart], z0m, obuk[ij]);
      }
  }

  return 0;
}

int cboundary_surface::surfm(double * restrict ustar, double * restrict obuk, 
                             double * restrict u, double * restrict ubot, double * restrict ugradbot, double * restrict ufluxbot, 
                             double * restrict v, double * restrict vbot, double * restrict vgradbot, double * restrict vfluxbot, 
                             double zsl, int bcbot)
{
  int ij,ijk,ii,jj,kk,kstart;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  // the surface value is known, calculate the flux and gradient
  if(bcbot == 0)
  {
    // first calculate the surface value
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        // interpolate the whole stability function rather than ustar or obuk
        ufluxbot[ij] = -(u[ijk]-ubot[ij])*0.5*(ustar[ij-ii]*fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*fm(zsl, z0m, obuk[ij]));
        vfluxbot[ij] = -(v[ijk]-vbot[ij])*0.5*(ustar[ij-jj]*fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*fm(zsl, z0m, obuk[ij]));
      }

    grid->boundary_cyclic2d(ufluxbot);
    grid->boundary_cyclic2d(vfluxbot);
  }
  // the flux is known, calculate the surface value and gradient
  else if(bcbot == 2)
  {
    // first redistribute ustar over the two flux components
    double u2, v2, vonu2,uonv2,ustaronu4,ustaronv4;
    const double minval = 1.e-2;

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        // minimize the wind at 0.01, thus the wind speed squared at 0.0001
        vonu2 = std::max(minval, 0.25*( std::pow(v[ijk-ii]-vbot[ij-ii], 2.) + std::pow(v[ijk-ii+jj]-vbot[ij-ii+jj], 2.)
                                      + std::pow(v[ijk   ]-vbot[ij   ], 2.) + std::pow(v[ijk   +jj]-vbot[ij   +jj], 2.)) );
        uonv2 = std::max(minval, 0.25*( std::pow(u[ijk-jj]-ubot[ij-jj], 2.) + std::pow(u[ijk+ii-jj]-ubot[ij+ii-jj], 2.)
                                      + std::pow(u[ijk   ]-ubot[ij   ], 2.) + std::pow(u[ijk+ii   ]-ubot[ij+ii   ], 2.)) );
        u2 = std::max(minval, std::pow(u[ijk]-ubot[ij], 2.) );
        v2 = std::max(minval, std::pow(v[ijk]-vbot[ij], 2.) );
        ustaronu4 = 0.5*(std::pow(ustar[ij-ii], 4.) + std::pow(ustar[ij], 4.));
        ustaronv4 = 0.5*(std::pow(ustar[ij-jj], 4.) + std::pow(ustar[ij], 4.));
        ufluxbot[ij] = -sign(u[ijk]) * std::pow(ustaronu4 / (1. + vonu2 / u2), 0.5);
        vfluxbot[ij] = -sign(v[ijk]) * std::pow(ustaronv4 / (1. + uonv2 / v2), 0.5);
      }

    grid->boundary_cyclic2d(ufluxbot);
    grid->boundary_cyclic2d(vfluxbot);
    
    // calculate the surface values
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        // interpolate the whole stability function rather than ustar or obuk
        ubot[ij] = ufluxbot[ij] / (0.5*(ustar[ij-ii]*fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + u[ijk];
        vbot[ij] = vfluxbot[ij] / (0.5*(ustar[ij-jj]*fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + v[ijk];
      }

    grid->boundary_cyclic2d(ubot);
    grid->boundary_cyclic2d(vbot);
  }

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      // use the linearly interpolated grad, rather than the MO grad,
      // to prevent giving unresolvable gradients to advection schemes
      // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0m*ustar[ij]) * phih(zsl/obuk[ij]);
      ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
      vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
    }

  return 0;
}

int cboundary_surface::surfs(double * restrict ustar, double * restrict obuk, double * restrict var,
                             double * restrict varbot, double * restrict vargradbot, double * restrict varfluxbot, 
                             double zsl, int bcbot)
{
  int ij,ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  // the surface value is known, calculate the flux and gradient
  if(bcbot == 0)
  {
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*fh(zsl, z0h, obuk[ij]);
        // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
        // use the linearly interpolated grad, rather than the MO grad,
        // to prevent giving unresolvable gradients to advection schemes
        vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
      }
  }
  else if(bcbot == 2)
  {
    // the flux is known, calculate the surface value and gradient
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        varbot[ij]     =  varfluxbot[ij] / (ustar[ij]*fh(zsl, z0h, obuk[ij])) + var[ijk];
        // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
        // use the linearly interpolated grad, rather than the MO grad,
        // to prevent giving unresolvable gradients to advection schemes
        vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
      }
  }

  return 0;
}

double cboundary_surface::calcobuk(double L, double du, double bfluxbot, double zsl)
{
  double L0;
  double Lstart, Lend;
  double fx, fxdif;

  int m = 0;
  int nlim = 10;

  const double Lmax = 1.e20;

  // avoid bfluxbot to be zero
  bfluxbot = std::max(dsmall, bfluxbot);

  // allow for one restart
  while(m <= 1)
  {
    // if L and bfluxbot are of the same sign, or the last calculation did not converge,
    // the stability has changed and the procedure needs to be reset
    if(L*bfluxbot >= 0.)
    {
      nlim = 200;
      if(bfluxbot >= 0.)
        L = -dsmall;
      else
        L = dsmall;
    }

    if(bfluxbot >= 0.)
      L0 = -dhuge;
    else
      L0 = dhuge;

    // TODO replace by value from buoyancy
    double gravitybeta = 9.81/300.;
    int n = 0;

    // exit on convergence or on iteration count
    while(std::abs((L - L0)/L0) > 0.001 && n < nlim && std::abs(L) < Lmax)
    {
      L0     = L;
      // fx     = Rib - zsl/L * (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L)) / std::pow(std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L), 2.);
      fx     = zsl/L + kappa*gravitybeta*zsl*bfluxbot / std::pow(du * fm(zsl, z0m, L), 3.);
      Lstart = L - 0.001*L;
      Lend   = L + 0.001*L;
      fxdif  = ( (zsl/Lend + kappa*gravitybeta*zsl*bfluxbot / std::pow(du * fm(zsl, z0m, Lend), 3.))
               - (zsl/Lstart + kappa*gravitybeta*zsl*bfluxbot / std::pow(du * fm(zsl, z0m, Lstart), 3.)) )
             / (Lend - Lstart);
      L      = L - fx/fxdif;
      ++n;
    }

    // convergence has been reached
    if(n < nlim && std::abs(L) < Lmax)
      break;
    // convergence has not been reached, procedure restarted once
    else
    {
      L = dsmall;
      ++m;
      nlim = 200;
    }
  }

  if(m > 1)
    std::printf("ERROR convergence has not been reached in Obukhov length calculation\n");

  return L;
}

inline double cboundary_surface::fm(double zsl, double z0m, double L)
{
  double fm;
  fm = kappa / (std::log(zsl/z0m) - psim(zsl/L) + psih(z0m/L));
  return fm;
}

inline double cboundary_surface::fh(double zsl, double z0h, double L)
{
  double fh;
  fh = kappa / (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L));
  return fh;
}

inline double cboundary_surface::psim(double zeta)
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

inline double cboundary_surface::psih(double zeta)
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

inline double cboundary_surface::phim(double zeta)
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

inline double cboundary_surface::phih(double zeta)
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
