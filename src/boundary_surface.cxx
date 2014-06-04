/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>    // std::count
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_surface.h"
#include "defines.h"
#include "thermo.h"
#include "model.h"
#include "master.h"
#include "cross.h"

#define NO_VELOCITY 0.
#define NO_OFFSET 0.

#define BC_DIRICHLET 0
#define BC_NEUMANN 1
#define BC_FLUX 2
#define BC_USTAR 3

// a sign function
inline double sign(double n) { return n > 0 ? 1 : (n < 0 ? -1 : 0);}

cboundary_surface::cboundary_surface(cmodel *modelin) : cboundary(modelin)
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
  int nerror = 0;

  nerror += processbcs(inputin);

  nerror += inputin->getItem(&z0m, "boundary", "z0m", "");
  nerror += inputin->getItem(&z0h, "boundary", "z0h", "");

  // Read list of cross sections
  nerror += inputin->getList(&crosslist , "boundary", "crosslist" , "");

  // copy all the boundary options and set the model ones to flux type
  surfmbcbot = mbcbot;
  mbcbot = BC_FLUX;

  // crash in case fixed gradient is prescribed
  if(surfmbcbot == BC_NEUMANN)
  {
    if(master->mpiid == 0) std::printf("ERROR neumann bc is not supported in surface model\n");
    ++nerror;
  }
  // read the ustar value only if fixed fluxes are prescribed
  else if(surfmbcbot == BC_USTAR)
    nerror += inputin->getItem(&ustarin, "boundary", "ustar", "");

  // process the scalars
  for(bcmap::const_iterator it=sbc.begin(); it!=sbc.end(); ++it)
  {
    surfsbcbot[it->first] = it->second->bcbot;
    //it->second->bcbot = BC_FLUX;

    // crash in case fixed gradient is prescribed
    if(surfsbcbot[it->first] == BC_NEUMANN)
    {
      if(master->mpiid == 0) std::printf("ERROR fixed gradient bc is not supported in surface model\n");
      ++nerror;
    }

    // crash in case of fixed momentum flux and dirichlet bc for scalar
    if(surfsbcbot[it->first] == BC_DIRICHLET && surfmbcbot == BC_USTAR)
    {
      if(master->mpiid == 0) std::printf("ERROR fixed ustar bc in combination with dirichlet bc for scalars is not supported\n");
      ++nerror;
    }
  }

  // check whether the prognostic thermo vars are of the same type
  std::vector<std::string> thermolist;
  model->thermo->getprogvars(&thermolist);

  std::vector<std::string>::const_iterator it = thermolist.begin();
  thermobc = surfsbcbot[*it];
  while(it != thermolist.end())
  {
    if(surfsbcbot[*it] != thermobc)
    {
      ++nerror;
      if(master->mpiid == 0) std::printf("ERROR all thermo variables need to have the same bc type\n");
    }
    ++it;
  }

  return nerror;
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

  // Cross sections
  allowedcrossvars.push_back("ustar");
  allowedcrossvars.push_back("obuk");

  // Check input list of cross variables (crosslist)
  std::vector<std::string>::iterator it=crosslist.begin();
  while(it != crosslist.end())
  {
    if(!std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it))
    {
      if(master->mpiid == 0) std::printf("WARNING field %s in [boundary][crosslist] is illegal\n", it->c_str());
      it = crosslist.erase(it);  // erase() returns iterator of next element..
    }
    else
      ++it;
  }

  return 0;
}

int cboundary_surface::execcross()
{
  int nerror = 0;

  for(std::vector<std::string>::iterator it=crosslist.begin(); it<crosslist.end(); ++it)
  {
    if(*it == "ustar")
      nerror += model->cross->crossplane(ustar, fields->s["tmp1"]->data, "ustar");
    else if(*it == "obuk")
      nerror += model->cross->crossplane(obuk,  fields->s["tmp1"]->data, "obuk");
  }  

  return nerror; 
}

int cboundary_surface::save(int iotime)
{
  char filename[256];

  std::sprintf(filename, "obuk.%07d", iotime);
  if(master->mpiid == 0) std::printf("Saving \"%s\" ... ", filename);
  if(grid->savexyslice(obuk, fields->s["tmp1"]->data, filename))
  {
    if(master->mpiid == 0) std::printf("FAILED\n");
    return 1;
  }
  else
    if(master->mpiid == 0) std::printf("OK\n");

  return 0;
}

int cboundary_surface::load(int iotime)
{
  char filename[256];

  std::sprintf(filename, "obuk.%07d", iotime);
  if(master->mpiid == 0) std::printf("Loading \"%s\" ... ", filename);
  if(grid->loadxyslice(obuk, fields->s["tmp1"]->data, filename))
  {
    if(master->mpiid == 0) std::printf("FAILED\n");
    return 1;
  }
  else
    if(master->mpiid == 0) std::printf("OK\n");

  grid->boundary_cyclic2d(obuk);

  return 0;
}

int cboundary_surface::setvalues()
{
  // grid transformation is properly taken into account by setting the databot and top values
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, surfmbcbot, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, surfmbcbot, NO_VELOCITY, fields->visc, grid->vtrans);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->vtrans);

  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setbc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, surfsbcbot[it->first], sbc[it->first]->bot, it->second->visc, NO_OFFSET);
    setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, NO_OFFSET);
  }

  // in case the momentum has a fixed ustar, set the value to that of the input
  if(surfmbcbot == BC_USTAR)
  {
    int ij,jj;
    jj = grid->icells;

    setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, 0, NO_VELOCITY, fields->visc, grid->utrans);
    setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, 0, NO_VELOCITY, fields->visc, grid->vtrans);

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
  if(model->thermo->getsw() == "0")
  {
    stability_neutral(ustar, obuk,
                      fields->u->data, fields->v->data,
                      fields->u->databot, fields->v->databot,
                      fields->sd["tmp1"]->data, grid->z);
  }
  else
  {
    // store the buoyancy in tmp1
    model->thermo->getbuoyancysurf(fields->sd["tmp1"]);
    stability(ustar, obuk, fields->sd["tmp1"]->datafluxbot,
              fields->u->data, fields->v->data, fields->sd["tmp1"]->data,
              fields->u->databot, fields->v->databot, fields->sd["tmp1"]->databot,
              fields->sd["tmp2"]->data, grid->z);
  }

  // calculate the surface value, gradient and flux depending on the chosen boundary condition
  surfm(ustar, obuk,
        fields->u->data, fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot,
        fields->v->data, fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot,
        grid->z[grid->kstart], surfmbcbot);

  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    surfs(ustar, obuk, it->second->data,
          it->second->databot, it->second->datagradbot, it->second->datafluxbot,
          grid->z[grid->kstart], surfsbcbot[it->first]);
  }

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
  double du2;
  //double utot, ubottot, du2;
  const double minval = 1.e-1;
  // first, interpolate the wind to the scalar location
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      // ubottot = std::pow(  0.5*(std::pow(ubot[ij], 2.) + std::pow(ubot[ij+ii], 2.))
      //                    + 0.5*(std::pow(vbot[ij], 2.) + std::pow(vbot[ij+jj], 2.)), 0.5);
      // utot    = std::pow(  0.5*(std::pow(u[ijk], 2.) + std::pow(u[ijk+ii], 2.))
      //                    + 0.5*(std::pow(v[ijk], 2.) + std::pow(v[ijk+jj], 2.)), 0.5);
      du2 = std::pow(0.5*(u[ijk] + u[ijk+ii]) - 0.5*(ubot[ij] + ubot[ij+ii]), 2)
          + std::pow(0.5*(v[ijk] + v[ijk+jj]) - 0.5*(vbot[ij] + vbot[ij+jj]), 2);
      // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
      // otherwise evisc at k = kstart blows up
      // dutot[ij] = std::max(std::abs(utot - ubottot), minval);
      dutot[ij] = std::max(std::pow(du2, 0.5), minval);
    }

  grid->boundary_cyclic2d(dutot);

  double db;

  // calculate Obukhov length
  // case 1: fixed buoyancy flux and fixed ustar
  if(surfmbcbot == BC_USTAR && thermobc == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        obuk[ij] = -std::pow(ustar[ij], 3.) / (kappa*bfluxbot[ij]);
      }
  }
  // case 2: fixed buoyancy surface value and free ustar
  else if(surfmbcbot == BC_DIRICHLET && thermobc == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        obuk [ij] = calcobuk_noslip_flux(obuk[ij], dutot[ij], bfluxbot[ij], z[kstart]);
        ustar[ij] = dutot[ij] * fm(z[kstart], z0m, obuk[ij]);
      }
  }
  else if(surfmbcbot == BC_DIRICHLET && thermobc == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        db = b[ijk] - bbot[ij];
        obuk [ij] = calcobuk_noslip_dirichlet(obuk[ij], dutot[ij], db, z[kstart]);
        ustar[ij] = dutot[ij] * fm(z[kstart], z0m, obuk[ij]);
      }
  }

  return 0;
}

int cboundary_surface::stability_neutral(double * restrict ustar, double * restrict obuk,
                                         double * restrict u    , double * restrict v   ,
                                         double * restrict ubot , double * restrict vbot,
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

  // set the Obukhov length to a very large negative number
  // case 1: fixed buoyancy flux and fixed ustar
  if(surfmbcbot == BC_USTAR && thermobc == BC_FLUX)
  {
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        obuk[ij] = -dbig;
      }
  }
  // case 2: fixed buoyancy surface value and free ustar
  else if(surfmbcbot == BC_DIRICHLET && thermobc == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        obuk [ij] = -dbig;
        ustar[ij] = dutot[ij] * fm(z[kstart], z0m, obuk[ij]);
      }
  }
  else if(surfmbcbot == BC_DIRICHLET && thermobc == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        obuk [ij] = -dbig;
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
  if(bcbot == BC_DIRICHLET)
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
  else if(bcbot == BC_USTAR)
  {
    // first redistribute ustar over the two flux components
    double u2,v2,vonu2,uonv2,ustaronu4,ustaronv4;
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
        ufluxbot[ij] = -sign(u[ijk]-ubot[ij]) * std::pow(ustaronu4 / (1. + vonu2 / u2), 0.5);
        vfluxbot[ij] = -sign(v[ijk]-vbot[ij]) * std::pow(ustaronv4 / (1. + uonv2 / v2), 0.5);
      }

    grid->boundary_cyclic2d(ufluxbot);
    grid->boundary_cyclic2d(vfluxbot);

    // CvH: I think that the problem is not closed, since both the fluxes and the surface values
    // of u and v are unknown. You have to assume a no slip in order to get the fluxes and therefore
    // should not update the surface values with those that belong to the flux. This procedure needs
    // to be checked more carefully.
    /*
    // calculate the surface values
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        // interpolate the whole stability function rather than ustar or obuk
        ubot[ij] = 0.;// ufluxbot[ij] / (0.5*(ustar[ij-ii]*fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + u[ijk];
        vbot[ij] = 0.;// vfluxbot[ij] / (0.5*(ustar[ij-jj]*fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + v[ijk];
      }

    grid->boundary_cyclic2d(ubot);
    grid->boundary_cyclic2d(vbot);
    */
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
  if(bcbot == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
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
  else if(bcbot == BC_FLUX)
  {
    // the flux is known, calculate the surface value and gradient
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        // if(ij=100) std::printf("CvH: ustar,fh, var[ijk]: %E, %E, %E\n", ustar[ij], obuk[ij], var[ijk]);
        varbot[ij] = varfluxbot[ij] / (ustar[ij]*fh(zsl, z0h, obuk[ij])) + var[ijk];
        // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
        // use the linearly interpolated grad, rather than the MO grad,
        // to prevent giving unresolvable gradients to advection schemes
        vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
      }
  }

  return 0;
}

double cboundary_surface::calcobuk_noslip_flux(double L, double du, double bfluxbot, double zsl)
{
  double L0;
  double Lstart, Lend;
  double fx, fxdif;

  int m = 0;
  int nlim = 10;

  const double Lmax = 1.e20;

  // avoid bfluxbot to be zero
  if(bfluxbot >= 0.)
    bfluxbot = std::max(dsmall, bfluxbot);
  else
    bfluxbot = std::min(-dsmall, bfluxbot);

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

    int n = 0;

    // exit on convergence or on iteration count
    while(std::abs((L - L0)/L0) > 0.001 && n < nlim && std::abs(L) < Lmax)
    {
      L0     = L;
      // fx     = Rib - zsl/L * (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L)) / std::pow(std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L), 2.);
      fx     = zsl/L + kappa*zsl*bfluxbot / std::pow(du * fm(zsl, z0m, L), 3.);
      Lstart = L - 0.001*L;
      Lend   = L + 0.001*L;
      fxdif  = ( (zsl/Lend + kappa*zsl*bfluxbot / std::pow(du * fm(zsl, z0m, Lend), 3.))
               - (zsl/Lstart + kappa*zsl*bfluxbot / std::pow(du * fm(zsl, z0m, Lstart), 3.)) )
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
double cboundary_surface::calcobuk_noslip_dirichlet(double L, double du, double db, double zsl)
{
  double L0;
  double Lstart, Lend;
  double fx, fxdif;

  int m = 0;
  int nlim = 10;

  const double Lmax = 1.e20;

  // avoid db to be zero
  if(db >= 0.)
    db = std::max(dsmall, db);
  else
    db = std::min(-dsmall, db);

  // allow for one restart
  while(m <= 1)
  {
    // if L and db are of different sign, or the last calculation did not converge,
    // the stability has changed and the procedure needs to be reset
    if(L*db <= 0.)
    {
      nlim = 200;
      if(db >= 0.)
        L = dsmall;
      else
        L = -dsmall;
    }

    if(db >= 0.)
      L0 = dhuge;
    else
      L0 = -dhuge;

    int n = 0;

    // exit on convergence or on iteration count
    while(std::abs((L - L0)/L0) > 0.001 && n < nlim && std::abs(L) < Lmax)
    {
      L0     = L;
      // fx     = Rib - zsl/L * (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L)) / std::pow(std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L), 2.);
      fx     = zsl/L - kappa*zsl*db*fh(zsl, z0h, L) / std::pow(du * fm(zsl, z0m, L), 2.);
      Lstart = L - 0.001*L;
      Lend   = L + 0.001*L;
      fxdif  = ( (zsl/Lend - kappa*zsl*db*fh(zsl, z0h, Lend) / std::pow(du * fm(zsl, z0m, Lend), 2.))
               - (zsl/Lstart - kappa*zsl*db*fh(zsl, z0h, Lstart) / std::pow(du * fm(zsl, z0m, Lstart), 2.)) )
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
    //exit(1);

  return L;
}

inline double cboundary_surface::fm(double zsl, double z0m, double L)
{
  double fm;
  fm = kappa / (std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L));
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
    psim = -2./3.*(zeta - 5./0.35) * std::exp(-0.35 * zeta) - zeta - (10./3.) / 0.35;
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
    psih  = (-2./3.) * (zeta-5./0.35) * std::exp(-0.35*zeta) - std::pow(1. + (2./3.) * zeta, 1.5) - (10./3.) / 0.35 + 1.;
  }
  return psih;
}

inline double cboundary_surface::phim(double zeta)
{
  double phim;
  if(zeta <= 0.)
  {
    // Businger-Dyer functions
    // phim  = (1. - 16. * zeta) ** (-0.25)
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
    // phih  = (1. - 16. * zeta) ** (-0.5)
    // Wilson functions
    phih = std::pow(1. + 7.9*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }
  else
    phih = 1. + 5.*zeta;

  return phih;
}
