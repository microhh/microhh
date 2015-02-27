/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "master.h"
#include "cross.h"
#include "most.h"

namespace
{
  // Size of the lookup table.
  const int nzL = 10000;
  inline double sign(double n) { return n > 0 ? 1 : (n < 0 ? -1 : 0); }
}

BoundarySurface::BoundarySurface(Model *modelin, Input *inputin) : Boundary(modelin, inputin)
{
  ustar = 0;
  obuk  = 0;
  nobuk = 0;

  #ifdef USECUDA
  ustar_g = 0;
  obuk_g  = 0;
  #endif
}

BoundarySurface::~BoundarySurface()
{
  delete[] ustar;
  delete[] obuk;
  delete[] nobuk;

  #ifdef USECUDA
  clearDevice();
  #endif
}

void BoundarySurface::create(Input *inputin)
{
  processTimeDep(inputin);

  // add variables to the statistics
  if(stats->getSwitch() == "1")
  {
    stats->addTimeSeries("ustar", "Surface friction velocity", "m s-1");
    stats->addTimeSeries("obuk", "Obukhov length", "m");
  }
}

void BoundarySurface::init(Input *inputin)
{
  // 1. Process the boundary conditions now all fields are registered
  processBcs(inputin);

  int nerror = 0;
  nerror += inputin->getItem(&z0m, "boundary", "z0m", "");
  nerror += inputin->getItem(&z0h, "boundary", "z0h", "");

  // Read list of cross sections
  nerror += inputin->getList(&crosslist , "boundary", "crosslist" , "");

  // copy all the boundary options and set the model ones to flux type
  // surfmbcbot = mbcbot;
  // mbcbot = FluxType;

  // crash in case fixed gradient is prescribed
  if(mbcbot == NeumannType)
  {
    master->printError("Neumann bc is not supported in surface model\n");
    ++nerror;
  }
  // read the ustar value only if fixed fluxes are prescribed
  else if(mbcbot == UstarType)
    nerror += inputin->getItem(&ustarin, "boundary", "ustar", "");

  // process the scalars
  for(BcMap::const_iterator it=sbc.begin(); it!=sbc.end(); ++it)
  {
    // surfsbcbot[it->first] = it->second->bcbot;
    // it->second->bcbot = FluxType;

    // crash in case fixed gradient is prescribed
    if(it->second->bcbot == NeumannType)
    {
      master->printError("fixed gradient bc is not supported in surface model\n");
      ++nerror;
    }

    // crash in case of fixed momentum flux and dirichlet bc for scalar
    if(it->second->bcbot == DirichletType && mbcbot == UstarType)
    {
      master->printError("fixed Ustar bc in combination with Dirichlet bc for scalars is not supported\n");
      ++nerror;
    }
  }

  // check whether the prognostic thermo vars are of the same type
  std::vector<std::string> thermolist;
  model->thermo->getProgVars(&thermolist);

  std::vector<std::string>::const_iterator it = thermolist.begin();

  // save the bc of the first thermo field in case thermo is enabled
  if(it != thermolist.end())
    thermobc = sbc[*it]->bcbot;

  while(it != thermolist.end())
  {
    if(sbc[*it]->bcbot != thermobc)
    {
      ++nerror;
      master->printError("all thermo variables need to have the same bc type\n");
    }
    ++it;
  }

  if(nerror)
    throw 1;

  // 2. Allocate the fields
  obuk  = new double[grid->ijcells];
  nobuk = new int   [grid->ijcells];
  ustar = new double[grid->ijcells];

  stats = model->stats;

  int ij,jj;
  jj = grid->icells;

  // initialize the obukhov length on a small number
  for(int j=0; j<grid->jcells; ++j)
    #pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij = i + j*jj;
      obuk[ij]  = constants::dsmall;
      nobuk[ij] = 0;
    }
 
  // Cross sections
  allowedcrossvars.push_back("ustar");
  allowedcrossvars.push_back("obuk");

  // Get global cross-list from cross.cxx
  std::vector<std::string> *crosslist_global = model->cross->getCrossList(); 

  // Check input list of cross variables (crosslist)
  std::vector<std::string>::iterator it2=crosslist_global->begin();
  while(it2 != crosslist_global->end())
  {
    if(std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it2))
    {
      // Remove variable from global list, put in local list
      crosslist.push_back(*it2);
      crosslist_global->erase(it2); // erase() returns iterator of next element..
    }
    else
      ++it2;
  }
}

void BoundarySurface::execCross()
{
  int nerror = 0;

  for(std::vector<std::string>::const_iterator it=crosslist.begin(); it<crosslist.end(); ++it)
  {
    if(*it == "ustar")
      nerror += model->cross->crossPlane(ustar, fields->atmp["tmp1"]->data, "ustar");
    else if(*it == "obuk")
      nerror += model->cross->crossPlane(obuk,  fields->atmp["tmp1"]->data, "obuk");
  }  

  if(nerror)
    throw 1;
}

void BoundarySurface::execStats(Mask *m)
{
  stats->calcMean2d(&m->tseries["obuk"].data , obuk , 0., fields->atmp["tmp4"]->databot, &stats->nmaskbot);
  stats->calcMean2d(&m->tseries["ustar"].data, ustar, 0., fields->atmp["tmp4"]->databot, &stats->nmaskbot);
}

void BoundarySurface::save(int iotime)
{
  char filename[256];

  std::sprintf(filename, "obuk.%07d", iotime);
  master->printMessage("Saving \"%s\" ... ", filename);
  if(grid->savexySlice(obuk, fields->atmp["tmp1"]->data, filename))
  {
    master->printMessage("FAILED\n");
    throw 1;
  }
  else
    master->printMessage("OK\n");
}

void BoundarySurface::load(int iotime)
{
  char filename[256];

  std::sprintf(filename, "obuk.%07d", iotime);
  master->printMessage("Loading \"%s\" ... ", filename);
  if(grid->loadxySlice(obuk, fields->atmp["tmp1"]->data, filename))
  {
    master->printMessage("FAILED\n");
    throw 1;
  }
  else
    master->printMessage("OK\n");

  grid->boundaryCyclic2d(obuk);
}

void BoundarySurface::setValues()
{
  const double noVelocity = 0.;
  const double noOffset = 0.;

  // grid transformation is properly taken into account by setting the databot and top values
  setBc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, noVelocity, fields->visc, grid->utrans);
  setBc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, noVelocity, fields->visc, grid->vtrans);

  setBc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, noVelocity, fields->visc, grid->utrans);
  setBc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, noVelocity, fields->visc, grid->vtrans);

  for(FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setBc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, noOffset);
    setBc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, noOffset);
  }

  // in case the momentum has a fixed ustar, set the value to that of the input
  if(mbcbot == UstarType)
  {
    int ij,jj;
    jj = grid->icells;

    setBc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, DirichletType, noVelocity, fields->visc, grid->utrans);
    setBc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, DirichletType, noVelocity, fields->visc, grid->vtrans);

    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
        {
          ij = i + j*jj;
          // limit ustar at 1e-4 to avoid zero divisions
          ustar[ij] = std::max(0.0001, ustarin);
        }
  }

  // Prepare the surface layer solver.
  zL_sl = new double[nzL];
  f_sl  = new double[nzL];

  const double dzL = 30./(nzL-1);
  zL_sl[0] = -20.;
  for (int n=1; n<nzL; ++n)
    zL_sl[n] = -20. + n*dzL;

  if(mbcbot == DirichletType && thermobc == FluxType)
  {
    const double zsl = grid->z[grid->kstart];
    for (int n=0; n<nzL; ++n)
      f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m, zsl/zL_sl[n]), 3);
  }
  else if(mbcbot == DirichletType && thermobc == DirichletType)
  {
    const double zsl = grid->z[grid->kstart];
    for (int n=0; n<nzL; ++n)
      f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m, zsl/zL_sl[n]), 2) / most::fh(zsl, z0h, zsl/zL_sl[n]);
  }
}

#ifndef USECUDA
void BoundarySurface::updateBcs()
{
  // start with retrieving the stability information
  if(model->thermo->getSwitch() == "0")
  {
    stabilityNeutral(ustar, obuk,
                     fields->u->data, fields->v->data,
                     fields->u->databot, fields->v->databot,
                     fields->atmp["tmp1"]->data, grid->z);
  }
  else
  {
    // store the buoyancy in tmp1
    model->thermo->getBuoyancySurf(fields->atmp["tmp1"]);
    stability(ustar, obuk, fields->atmp["tmp1"]->datafluxbot,
              fields->u->data,    fields->v->data,    fields->atmp["tmp1"]->data,
              fields->u->databot, fields->v->databot, fields->atmp["tmp1"]->databot,
              fields->atmp["tmp2"]->data, grid->z);
  }

  // calculate the surface value, gradient and flux depending on the chosen boundary condition
  surfm(ustar, obuk,
        fields->u->data, fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot,
        fields->v->data, fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot,
        grid->z[grid->kstart], mbcbot);

  for(FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    surfs(ustar, obuk, it->second->data,
          it->second->databot, it->second->datagradbot, it->second->datafluxbot,
          grid->z[grid->kstart], sbc[it->first]->bcbot);
  }
}
#endif

void BoundarySurface::stability(double * restrict ustar, double * restrict obuk, double * restrict bfluxbot,
                               double * restrict u    , double * restrict v   , double * restrict b       ,
                               double * restrict ubot , double * restrict vbot, double * restrict bbot    ,
                               double * restrict dutot, double * restrict z)
{
  int ij,ijk,ii,jj,kk,kstart;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

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
      // ubottot = std::pow(  0.5*(std::pow(ubot[ij], 2) + std::pow(ubot[ij+ii], 2))
      //                    + 0.5*(std::pow(vbot[ij], 2) + std::pow(vbot[ij+jj], 2)), 0.5);
      // utot    = std::pow(  0.5*(std::pow(u[ijk], 2) + std::pow(u[ijk+ii], 2))
      //                    + 0.5*(std::pow(v[ijk], 2) + std::pow(v[ijk+jj], 2)), 0.5);
      du2 = std::pow(0.5*(u[ijk] + u[ijk+ii]) - 0.5*(ubot[ij] + ubot[ij+ii]), 2)
          + std::pow(0.5*(v[ijk] + v[ijk+jj]) - 0.5*(vbot[ij] + vbot[ij+jj]), 2);
      // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
      // otherwise evisc at k = kstart blows up
      // dutot[ij] = std::max(std::abs(utot - ubottot), minval);
      dutot[ij] = std::max(std::pow(du2, 0.5), minval);
    }

  grid->boundaryCyclic2d(dutot);

  double db;

  // calculate Obukhov length
  // case 1: fixed buoyancy flux and fixed ustar
  if(mbcbot == UstarType && thermobc == FluxType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        obuk[ij] = -std::pow(ustar[ij], 3) / (constants::kappa*bfluxbot[ij]);
      }
  }
  // case 2: fixed buoyancy surface value and free ustar
  else if(mbcbot == DirichletType && thermobc == FluxType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        obuk [ij] = calcObukNoslipFlux(zL_sl, f_sl, nobuk[ij], dutot[ij], bfluxbot[ij], z[kstart]);
        ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m, obuk[ij]);
      }
  }
  else if(mbcbot == DirichletType && thermobc == DirichletType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        db = b[ijk] - bbot[ij];
        obuk [ij] = calcObukNoslipDirichlet(zL_sl, f_sl, nobuk[ij], dutot[ij], db, z[kstart]);
        ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m, obuk[ij]);
      }
  }
}

void BoundarySurface::stabilityNeutral(double * restrict ustar, double * restrict obuk,
                                       double * restrict u    , double * restrict v   ,
                                       double * restrict ubot , double * restrict vbot,
                                       double * restrict dutot, double * restrict z)
{
  int ij,ijk,ii,jj,kk,kstart;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  kstart = grid->kstart;

  // calculate total wind
  double du2;
  //double utot, ubottot;
  const double minval = 1.e-1;
  // first, interpolate the wind to the scalar location
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      //ubottot = std::pow(  0.5*(std::pow(ubot[ij], 2.) + std::pow(ubot[ij+ii], 2.))
      //                   + 0.5*(std::pow(vbot[ij], 2.) + std::pow(vbot[ij+jj], 2.)), 0.5);
      //utot    = std::pow(  0.5*(std::pow(u[ijk], 2.) + std::pow(u[ijk+ii], 2.))
      //                   + 0.5*(std::pow(v[ijk], 2.) + std::pow(v[ijk+jj], 2.)), 0.5);
      du2 = std::pow(0.5*(u[ijk] + u[ijk+ii]) - 0.5*(ubot[ij] + ubot[ij+ii]), 2)
          + std::pow(0.5*(v[ijk] + v[ijk+jj]) - 0.5*(vbot[ij] + vbot[ij+jj]), 2);
      // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
      // otherwise evisc at k = kstart blows up
      //dutot[ij] = std::max(std::abs(utot - ubottot), minval);
      dutot[ij] = std::max(std::pow(du2, 0.5), minval);
    }

  grid->boundaryCyclic2d(dutot);

  // set the Obukhov length to a very large negative number
  // case 1: fixed buoyancy flux and fixed ustar
  if(mbcbot == UstarType && thermobc == FluxType)
  {
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        obuk[ij] = -constants::dbig;
      }
  }
  // case 2: fixed buoyancy surface value and free ustar
  else if(mbcbot == DirichletType && thermobc == FluxType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        obuk [ij] = -constants::dbig;
        ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m, obuk[ij]);
      }
  }
  else if(mbcbot == DirichletType && thermobc == DirichletType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        obuk [ij] = -constants::dbig;
        ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m, obuk[ij]);
      }
  }
}

void BoundarySurface::surfm(double * restrict ustar, double * restrict obuk, 
                             double * restrict u, double * restrict ubot, double * restrict ugradbot, double * restrict ufluxbot, 
                             double * restrict v, double * restrict vbot, double * restrict vgradbot, double * restrict vfluxbot, 
                             double zsl, int bcbot)
{
  int ij,ijk,ii,jj,kk,kstart;

  ii = 1;
  jj = grid->icells;
  kk = grid->ijcells;

  kstart = grid->kstart;

  // the surface value is known, calculate the flux and gradient
  if(bcbot == DirichletType)
  {
    // first calculate the surface value
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        // interpolate the whole stability function rather than ustar or obuk
        ufluxbot[ij] = -(u[ijk]-ubot[ij])*0.5*(ustar[ij-ii]*most::fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*most::fm(zsl, z0m, obuk[ij]));
        vfluxbot[ij] = -(v[ijk]-vbot[ij])*0.5*(ustar[ij-jj]*most::fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*most::fm(zsl, z0m, obuk[ij]));
      }

    grid->boundaryCyclic2d(ufluxbot);
    grid->boundaryCyclic2d(vfluxbot);
  }
  // the flux is known, calculate the surface value and gradient
  else if(bcbot == UstarType)
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
        vonu2 = std::max(minval, 0.25*( std::pow(v[ijk-ii]-vbot[ij-ii], 2) + std::pow(v[ijk-ii+jj]-vbot[ij-ii+jj], 2)
                                      + std::pow(v[ijk   ]-vbot[ij   ], 2) + std::pow(v[ijk   +jj]-vbot[ij   +jj], 2)) );
        uonv2 = std::max(minval, 0.25*( std::pow(u[ijk-jj]-ubot[ij-jj], 2) + std::pow(u[ijk+ii-jj]-ubot[ij+ii-jj], 2)
                                      + std::pow(u[ijk   ]-ubot[ij   ], 2) + std::pow(u[ijk+ii   ]-ubot[ij+ii   ], 2)) );
        u2 = std::max(minval, std::pow(u[ijk]-ubot[ij], 2) );
        v2 = std::max(minval, std::pow(v[ijk]-vbot[ij], 2) );
        ustaronu4 = 0.5*(std::pow(ustar[ij-ii], 4) + std::pow(ustar[ij], 4));
        ustaronv4 = 0.5*(std::pow(ustar[ij-jj], 4) + std::pow(ustar[ij], 4));
        ufluxbot[ij] = -sign(u[ijk]-ubot[ij]) * std::pow(ustaronu4 / (1. + vonu2 / u2), 0.5);
        vfluxbot[ij] = -sign(v[ijk]-vbot[ij]) * std::pow(ustaronv4 / (1. + uonv2 / v2), 0.5);
      }

    grid->boundaryCyclic2d(ufluxbot);
    grid->boundaryCyclic2d(vfluxbot);

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

    grid->boundaryCyclic2d(ubot);
    grid->boundaryCyclic2d(vbot);
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
}

void BoundarySurface::surfs(double * restrict ustar, double * restrict obuk, double * restrict var,
                            double * restrict varbot, double * restrict vargradbot, double * restrict varfluxbot, 
                            double zsl, int bcbot)
{
  int ij,ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->ijcells;

  kstart = grid->kstart;

  // the surface value is known, calculate the flux and gradient
  if(bcbot == DirichletType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h, obuk[ij]);
        // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
        // use the linearly interpolated grad, rather than the MO grad,
        // to prevent giving unresolvable gradients to advection schemes
        vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
      }
  }
  else if(bcbot == FluxType)
  {
    // the flux is known, calculate the surface value and gradient
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        // if(ij=100) std::printf("CvH: ustar,fh, var[ijk]: %E, %E, %E\n", ustar[ij], obuk[ij], var[ijk]);
        varbot[ij] = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h, obuk[ij])) + var[ijk];
        // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
        // use the linearly interpolated grad, rather than the MO grad,
        // to prevent giving unresolvable gradients to advection schemes
        vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
      }
  }
}

double BoundarySurface::calcObukNoslipFlux(const double* const restrict zL, const double* const restrict f,
                                           int& n,
                                           const double du, const double bfluxbot, const double zsl)
{
  // Calculate the appropriate Richardson number.
  const double Ri = -constants::kappa * bfluxbot * zsl / std::pow(du, 3);

  // Determine search direction.
  if ( (f[n]-Ri) > 0)
  {
    while (f[n]-Ri > 0)
    {
      if ( (f[n]-Ri) < 0 || n == 0 )
        break;
      else
        --n;
    }
  }
  else
  {
    while ( (f[n]-Ri) < 0)
    {
      ++n;
      if ( (f[n]-Ri) > 0 || n == nzL-1 )
        break;
    }
  }

  double zL0;
  if (n == 0)
  {
    zL0 = zL[n];
    master->printWarning("z/L range too limited on unstable side\n");
  }
  else if (n == nzL)
  {
    zL0 = zL[n];
    master->printWarning("z/L range too limited on stable side\n");
  }
  else
  {
    // Linearly interpolate to the correct value of z/L.
    zL0 = zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);
    master->printMessage("%E, %E, %E\n", zL, f[n-1], f[n]);
  }

  return zsl/zL0;
}

double BoundarySurface::calcObukNoslipDirichlet(const double* const restrict zL, const double* const restrict f,
                                                int& n,
                                                const double du, const double db, const double zsl)
{
  // Calculate the appropriate Richardson number.
  const double Ri = constants::kappa * db * zsl / std::pow(du, 2);

  // Determine search direction.
  if ( (f[n]-Ri) > 0)
  {
    while (f[n-1]-Ri > 0 && n > 0)
    {
      --n;
    }
  }
  else
  {
    while ( (f[n]-Ri) < 0 && n < nzL-1)
    {
      ++n;
    }
  }

  double zL0;
  if (n == 0)
  {
    zL0 = zL[n];
    master->printWarning("z/L range too limited on unstable side\n");
  }
  else if (n == nzL)
  {
    zL0 = zL[n];
    master->printWarning("z/L range too limited on stable side\n");
  }
  else
  {
    // Linearly interpolate to the correct value of z/L.
    zL0 = zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);
    // master->printMessage("%d, %E, %E\n", n, f[n-1]-Ri, f[n]-Ri);
  }

  // master->printMessage("%E, %E\n", zL0, zsl/zL0);
  return zsl/zL0;
}
