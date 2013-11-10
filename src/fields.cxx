#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "defines.h"

#define NO_OFFSET 0.

cfields::cfields(cgrid *gridin, cmpi *mpiin)
{
  grid = gridin;
  mpi  = mpiin;

  allocated = false;
}

cfields::~cfields()
{
  if(allocated)
  {
    // DEALLOCATE ALL THE FIELDS
    // deallocate the prognostic velocity fields
    for(fieldmap::iterator it=mp.begin(); it!=mp.end(); ++it)
      delete it->second;

    // deallocate the velocity tendency fields
    for(fieldmap::iterator it=mt.begin(); it!=mt.end(); ++it)
      delete it->second;

    // deallocate the prognostic scalar fields
    for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
      delete it->second;

    // deallocate the scalar tendency fields
    for(fieldmap::iterator it=st.begin(); it!=st.end(); ++it)
      delete it->second;

    // deallocate the diagnostic scalars
    for(fieldmap::iterator it=sd.begin(); it!=sd.end(); ++it)
      delete it->second;
  }
}

int cfields::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&visc , "fields", "visc", "" );

  // LES
  n += inputin->getItem(&tPr, "fields", "tPr", "", 1./3.);

  // read the name of the passive scalars
  std::vector<std::string> slist;
  n += inputin->getList(&slist, "fields", "slist", "");

  // initialize the scalars
  for(std::vector<std::string>::iterator it = slist.begin(); it!=slist.end(); ++it)
  {
    if(initpfld(*it))
      return 1;
    n += inputin->getItem(&sp[*it]->visc, "fields", "svisc", *it);
  }

  // initialize the basic set of fields
  n += initmomfld(u, ut, "u");
  n += initmomfld(v, vt, "v");
  n += initmomfld(w, wt, "w");
  n += initdfld("p");
  n += initdfld("tmp1");
  n += initdfld("tmp2");

  if(n > 0)
    return 1;

  return 0;
}

int cfields::init()
{
  if(mpi->mpiid == 0) std::printf("Initializing fields\n");

  int n = 0;

  // ALLOCATE ALL THE FIELDS
  // allocate the prognostic velocity fields
  for(fieldmap::iterator it=mp.begin(); it!=mp.end(); ++it)
    n += it->second->init();

  // allocate the velocity tendency fields
  for(fieldmap::iterator it=mt.begin(); it!=mt.end(); ++it)
    n += it->second->init();

  // allocate the prognostic scalar fields
  for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
    n += it->second->init();

  // allocate the scalar tendency fields
  for(fieldmap::iterator it=st.begin(); it!=st.end(); ++it)
    n += it->second->init();

  // allocate the diagnostic scalars
  for(fieldmap::iterator it=sd.begin(); it!=sd.end(); ++it)
    n += it->second->init();

  if(n > 0)
    return 1;

  allocated = true;

  return 0;
}

int cfields::initmomfld(cfield3d *&fld, cfield3d *&fldt, std::string fldname)
{
  if (mp.find(fldname)!=mp.end())
  {
    std::printf("ERROR \"%s\" already exists\n", fldname.c_str());
    return 1;
  }

  // add a new prognostic momentum variable
  mp[fldname] = new cfield3d(grid, mpi, fldname);

  // add a new tendency for momentum variable
  std::string fldtname = fldname + "t";
  mt[fldname] = new cfield3d(grid, mpi, fldtname);

  // TODO remove these from the model?
  fld  = mp[fldname];
  fldt = mt[fldname];

  // add the prognostic variable and its tendency to the collection
  // of all fields and tendencies
  ap[fldname] = mp[fldname];
  at[fldname] = mt[fldname];

  return 0;
}

int cfields::initpfld(std::string fldname)
{
  if (s.find(fldname)!=s.end())
  {
    std::printf("ERROR \"%s\" already exists\n", fldname.c_str());
    return 1;
  }
  
  // add a new scalar variable
  sp[fldname] = new cfield3d(grid, mpi, fldname);

  // add a new tendency for scalar variable
  std::string fldtname = fldname + "t";
  st[fldname] = new cfield3d(grid, mpi, fldtname);

  // add the prognostic variable and its tendency to the collection
  // of all fields and tendencies
  s [fldname] = sp[fldname];
  ap[fldname] = sp[fldname];
  at[fldname] = st[fldname];
  
  return 0;
}

int cfields::initdfld(std::string fldname)
{
  if (s.find(fldname)!=s.end())
  {
    std::printf("ERROR \"%s\" already exists\n", fldname.c_str());
    return 1;
  }

  sd[fldname] = new cfield3d(grid, mpi, fldname );
  s [fldname] = sd[fldname];

  return 0;  
}

int cfields::create(cinput *inputin)
{
  if(mpi->mpiid == 0) std::printf("Creating fields\n");
  
  int n = 0;
  
  // Randomnize the momentum
  for(fieldmap::iterator it=mp.begin(); it!=mp.end(); ++it)
    n += randomnize(inputin, it->first, it->second->data);
  
  // Randomnize the scalars
  for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
    n += randomnize(inputin, it->first, it->second->data);
  
  // Add Vortices
  n += addvortexpair(inputin);
  
  // Add the mean profiles to the fields
  n += addmeanprofile(inputin, "u", mp["u"]->data, grid->u);
  n += addmeanprofile(inputin, "v", mp["v"]->data, grid->v);
 
  for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
    n += addmeanprofile(inputin, it->first, it->second->data, 0.);
  
  // set w equal to zero at the boundaries, just to be sure
  int lbot = grid->kstart*grid->icells*grid->jcells;
  int ltop = grid->kend  *grid->icells*grid->jcells;
  for(int l=0; l<grid->icells*grid->jcells; ++l)
  {
    w->data[lbot+l] = 0.;
    w->data[ltop+l] = 0.;
  }
  
  return (n>0);
}

int cfields::randomnize(cinput *inputin, std::string fld, double * restrict data)
{
  int n = 0;

  // set mpiid as random seed to avoid having the same field at all procs
  int static seed = 0;

  if(!seed)
  {
    n += inputin->getItem(&seed, "fields", "rndseed" , "", 2);
    seed += mpi->mpiid;
    std::srand(seed);
  }
  
  int ijk,jj,kk;
  int kendrnd;
  double rndfac;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  // look up the specific randomnizer variables
  n += inputin->getItem(&rndamp , "fields", "rndamp" , fld, 0.);
  n += inputin->getItem(&rndz   , "fields", "rndz"   , fld, 0.);
  n += inputin->getItem(&rndbeta, "fields", "rndbeta", fld, 0.);

  // find the location of the randomizer height
  kendrnd = grid->kstart;
  while(grid->zh[kendrnd+1] < rndz)
  {
    ++kendrnd;
  }

  if(kendrnd > grid->kend)
  {
    printf("ERROR: Randomnizer height rndz (%f) higher than domain top (%f)\n", grid->z[kendrnd],grid->zsize);
    return 1;
  }
  
  if((int) kendrnd == grid->kstart)
  {
    kendrnd = grid->kend;
  }

  for(int k=grid->kstart; k<kendrnd; ++k)
  {
    rndfac = std::pow((rndz-grid->z [k])/rndz, rndbeta);
    for(int j=grid->jstart; j<grid->jend; ++j)
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] = rndfac  * rndamp * ((double) std::rand() / (double) RAND_MAX - 0.5);
      }
  }

  return (n>0);
}

int cfields::addvortexpair(cinput *inputin)
{
  // add a double vortex to the initial conditions
  const double pi = std::acos((double)-1.);
  int ijk, jj, kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  int n = 0;
    
  // optional parameters
  n += inputin->getItem(&nvortexpair, "fields", "nvortexpair", "", 0    );
  n += inputin->getItem(&vortexamp  , "fields", "vortexamp"  , "", 1.e-3);
  n += inputin->getItem(&vortexaxis , "fields", "vortexaxis" , "", 1    );

  if(nvortexpair > 0)
  {
    if(vortexaxis == 0)
      for(int k=grid->kstart; k<grid->kend; ++k)
        for(int j=grid->jstart; j<grid->jend; ++j)
          for(int i=grid->istart; i<grid->iend; ++i)
          {
            ijk = i + j*jj + k*kk;
            u->data[ijk] +=  vortexamp*std::sin(nvortexpair*2.*pi*(grid->xh[i])/grid->xsize)*std::cos(pi*grid->z [k]/grid->zsize);
            w->data[ijk] += -vortexamp*std::cos(nvortexpair*2.*pi*(grid->x [i])/grid->xsize)*std::sin(pi*grid->zh[k]/grid->zsize);
          }
    else if(vortexaxis == 1)
      for(int k=grid->kstart; k<grid->kend; ++k)
        for(int j=grid->jstart; j<grid->jend; ++j)
          for(int i=grid->istart; i<grid->iend; ++i)
          {
            ijk = i + j*jj + k*kk;
            v->data[ijk] +=  vortexamp*std::sin(nvortexpair*2.*pi*(grid->yh[j])/grid->ysize)*std::cos(pi*grid->z [k]/grid->zsize);
            w->data[ijk] += -vortexamp*std::cos(nvortexpair*2.*pi*(grid->y [j])/grid->ysize)*std::sin(pi*grid->zh[k]/grid->zsize);
          }
  }
  
  return (n>0);
}

int cfields::addmeanprofile(cinput *inputin, std::string fld, double * restrict data, double offset)
{
  int ijk, jj, kk;
  double proftemp[grid->kmax];

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  if(inputin->getProf(proftemp, fld, grid->kmax))
    return 1;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] += proftemp[k-grid->kstart] - offset;
      }
      
  return 0;
}

int cfields::load(int n)
{
  int nerror = 0;

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    // the offset is kept at zero, otherwise bitwise identical restarts is not possible
    char filename[256];
    std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
    if(mpi->mpiid == 0) std::printf("Loading \"%s\" ... ", filename);
    if(grid->loadfield3d(it->second->data, sd["tmp1"]->data, sd["tmp2"]->data, filename, NO_OFFSET))
    {
      if(mpi->mpiid == 0) std::printf("FAILED\n");
      ++nerror;
    }
    else
    {
      if(mpi->mpiid == 0) std::printf("OK\n");
    }  
  }

  return nerror;
}

int cfields::save(int n)
{
  int nerror = 0;
  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", it->second->name.c_str(), n);
    if(mpi->mpiid == 0) std::printf("Saving \"%s\" ... ", filename);

    // the offset is kept at zero, because otherwise bitwise identical restarts is not possible
    if(grid->savefield3d(it->second->data, sd["tmp1"]->data, sd["tmp2"]->data, filename, NO_OFFSET))
    {
      if(mpi->mpiid == 0) std::printf("FAILED\n");
      ++nerror;
    }  
    else
    {
      if(mpi->mpiid == 0) std::printf("OK\n");
    }
  }

  return nerror;
}

double cfields::checkmom()
{
  return calcmom_2nd(u->data, v->data, w->data, grid->dz);
}

double cfields::checktke()
{
  return calctke_2nd(u->data, v->data, w->data, grid->dz);
}

double cfields::checkmass()
{
  // CvH for now, do the mass check on the first scalar... Do we want to change this?
  fieldmap::iterator itProg=sp.begin();
  if(sp.begin() != sp.end())
    return calcmass(itProg->second->data, grid->dz);
  else
    return 0.;
}

double cfields::calcmass(double * restrict s, double * restrict dz)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double mass = 0;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        mass += s[ijk]*dz[k];
      }

  grid->getsum(&mass);

  mass /= (grid->itot*grid->jtot*grid->zsize);

  return mass;
}

double cfields::calcmom_2nd(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double momentum;
  momentum = 0;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        momentum += (interp2(u[ijk], u[ijk+ii]) + interp2(v[ijk], v[ijk+jj]) + interp2(w[ijk], w[ijk+kk]))*dz[k];
      }

  grid->getsum(&momentum);

  momentum /= (grid->itot*grid->jtot*grid->zsize);

  return momentum;
}

double cfields::calctke_2nd(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double tke = 0;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        tke += ( interp2(u[ijk]*u[ijk], u[ijk+ii]*u[ijk+ii]) 
               + interp2(v[ijk]*v[ijk], v[ijk+jj]*v[ijk+jj]) 
               + interp2(w[ijk]*w[ijk], w[ijk+kk]*w[ijk+kk]))*dz[k];
      }

  grid->getsum(&tke);

  tke /= (grid->itot*grid->jtot*grid->zsize);
  tke *= 0.5;

  return tke;
}

inline double cfields::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

