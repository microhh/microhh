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

#ifndef PARALLEL
#include <fftw3.h>
#include <cstdio>
#include "master.h"
#include "grid.h"
#include "defines.h"

// MPI functions
int cgrid::initmpi()
{
  mpitypes = true;
  return 0;
} 

int cgrid::exitmpi()
{
  return 0;
}

int cgrid::boundary_cyclic(double * restrict data)
{
  int ncount = 1;
  int ijk0,ijk1,jj,kk;

  jj = icells;
  kk = icells*jcells;

  // first, east west boundaries
  for(int k=0; k<kcells; k++)
    for(int j=0; j<jcells; j++)
#pragma ivdep
      for(int i=0; i<igc; i++)
      {
        ijk0 = i          + j*jj + k*kk;
        ijk1 = iend-igc+i + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<kcells; k++)
    for(int j=0; j<jcells; j++)
#pragma ivdep
      for(int i=0; i<igc; i++)
      {
        ijk0 = i+iend   + j*jj + k*kk;
        ijk1 = i+istart + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  // if the run is 3D, apply the BCs
  if(jtot > 1)
  {
    // second, send and receive the ghost cells in the north-south direction
    for(int k=0; k<kcells; k++)
      for(int j=0; j<jgc; j++)
#pragma ivdep
        for(int i=0; i<icells; i++)
        {
          ijk0 = i + j           *jj + k*kk;
          ijk1 = i + (jend-jgc+j)*jj + k*kk;
          data[ijk0] = data[ijk1];
        }

    for(int k=0; k<kcells; k++)
      for(int j=0; j<jgc; j++)
#pragma ivdep
        for(int i=0; i<icells; i++)
        {
          ijk0 = i + (j+jend  )*jj + k*kk;
          ijk1 = i + (j+jstart)*jj + k*kk;
          data[ijk0] = data[ijk1];
        }
  }
  // in case of 2D, fill all the ghost cells with the current value
  else
  {
    // 2d essential variables
    int ijkref,ijknorth,ijksouth,jj,kk;

    jj = icells;
    kk = icells*jcells;

    for(int k=kstart; k<kend; k++)
      for(int j=0; j<jgc; j++)
#pragma ivdep
        for(int i=0; i<icells; i++)
        {
          ijkref   = i + jstart*jj   + k*kk;
          ijknorth = i + j*jj        + k*kk;
          ijksouth = i + (jend+j)*jj + k*kk;
          data[ijknorth] = data[ijkref];
          data[ijksouth] = data[ijkref];
        }
  }

  return 0;
}

int cgrid::boundary_cyclic2d(double * restrict data)
{
  int ij0,ij1,jj;

  jj = icells;

  // first, east west boundaries
  for(int j=0; j<jcells; j++)
#pragma ivdep
    for(int i=0; i<igc; i++)
    {
      ij0 = i          + j*jj;
      ij1 = iend-igc+i + j*jj;
      data[ij0] = data[ij1];
    }

  for(int j=0; j<jcells; j++)
#pragma ivdep
    for(int i=0; i<igc; i++)
    {
      ij0 = i+iend   + j*jj;
      ij1 = i+istart + j*jj;
      data[ij0] = data[ij1];
    }

  // if the run is 3D, apply the BCs
  if(jtot > 1)
  {
    // second, send and receive the ghost cells in the north-south direction
    for(int j=0; j<jgc; j++)
#pragma ivdep
      for(int i=0; i<icells; i++)
      {
        ij0 = i + j           *jj;
        ij1 = i + (jend-jgc+j)*jj;
        data[ij0] = data[ij1];
      }

    for(int j=0; j<jgc; j++)
#pragma ivdep
      for(int i=0; i<icells; i++)
      {
        ij0 = i + (j+jend  )*jj;
        ij1 = i + (j+jstart)*jj;
        data[ij0] = data[ij1];
      }
  }
  // in case of 2D, fill all the ghost cells with the current value
  else
  {
    // 2d essential variables
    int ijref,ijnorth,ijsouth,jj;

    jj = icells;

    for(int j=0; j<jgc; j++)
#pragma ivdep
      for(int i=0; i<icells; i++)
      {
        ijref   = i + jstart*jj;
        ijnorth = i + j*jj;
        ijsouth = i + (jend+j)*jj;
        data[ijnorth] = data[ijref];
        data[ijsouth] = data[ijref];
      }
  }

  return 0;
}

int cgrid::transposezx(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposexz(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposexy(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposeyx(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposeyz(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::transposezy(double * restrict ar, double * restrict as)
{
  int ijk,jj,kk;
  jj = imax;
  kk = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk = i + j*jj + k*kk;
        ar[ijk] = as[ijk];
      }

  return 0;
}

int cgrid::getmax(double *var)
{
  return 0;
}

int cgrid::getsum(double *var)
{
  return 0;
}

int cgrid::getprof(double *prof, int kcellsin)
{
  return 0;
}

// IO functions
int cgrid::save()
{
  // SAVE THE GRID
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  pFile = fopen(filename, "wb");
  std::printf("Saving \"%s\" ... ", filename);

  if(pFile == NULL)
  {
    std::printf("FAILED\n");
    return 1;
  }
  else
    std::printf("OK\n");

  fwrite(&x [istart], sizeof(double), itot, pFile);
  fwrite(&xh[istart], sizeof(double), itot, pFile);
  fwrite(&y [jstart], sizeof(double), jtot, pFile);
  fwrite(&yh[jstart], sizeof(double), jtot, pFile);
  fwrite(&z [kstart], sizeof(double), ktot, pFile);
  fwrite(&zh[kstart], sizeof(double), ktot, pFile);
  fclose(pFile);

  // SAVE THE FFTW PLAN IN ORDER TO ENSURE BITWISE IDENTICAL RESTARTS
  // use the FFTW3 many interface in order to reduce function call overhead
  int rank = 1;
  int ni[] = {itot};
  int nj[] = {jtot};
  int istride = 1;
  int jstride = iblock;
  int idist = itot;
  int jdist = 1;
  fftw_r2r_kind kindf[] = {FFTW_R2HC};
  fftw_r2r_kind kindb[] = {FFTW_HC2R};
  iplanf = fftw_plan_many_r2r(1, ni, jmax, fftini, ni, istride, idist,
                              fftouti, ni, istride, idist, kindf, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_many_r2r(1, ni, jmax, fftini, ni, istride, idist,
                              fftouti, ni, istride, idist, kindb, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_many_r2r(1, nj, iblock, fftinj, nj, jstride, jdist,
                              fftoutj, nj, jstride, jdist, kindf, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_many_r2r(1, nj, iblock, fftinj, nj, jstride, jdist,
                              fftoutj, nj, jstride, jdist, kindb, FFTW_EXHAUSTIVE);

  fftwplan = true;

  if(master->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", "fftwplan", 0);

    std::printf("Saving \"%s\" ... ", filename);

    int n = fftw_export_wisdom_to_filename(filename);
    if(n == 0)
    {
      std::printf("FAILED\n");
      return 1;
    }
    else
      std::printf("OK\n");
  }

  return 0;
}

int cgrid::load()
{
  // LOAD THE GRID
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  pFile = fopen(filename, "rb");
  std::printf("Loading \"%s\" ... ", filename);

  if(pFile == NULL)
  {
    std::printf("FAILED\n");
    return 1;
  }
  else
    std::printf("OK\n");

  fread(&x [istart], sizeof(double), itot, pFile);
  fread(&xh[istart], sizeof(double), itot, pFile);
  fread(&y [jstart], sizeof(double), jtot, pFile);
  fread(&yh[jstart], sizeof(double), jtot, pFile);
  fread(&z [kstart], sizeof(double), ktot, pFile);
  fread(&zh[kstart], sizeof(double), ktot, pFile);
  fclose(pFile);

  // calculate the missing coordinates
  calculate();

  // LOAD THE FFTW PLAN
  std::sprintf(filename, "%s.%07d", "fftwplan", 0);

  if(master->mpiid == 0)
    std::printf("Loading \"%s\" ... ", filename);

  int n = fftw_import_wisdom_from_filename(filename);
  if(n == 0)
  {
    std::printf("FAILED\n");
    return 1;
  }
  else
    std::printf("OK\n");

  // use the FFTW3 many interface in order to reduce function call overhead
  int rank = 1;
  int ni[] = {itot};
  int nj[] = {jtot};
  int istride = 1;
  int jstride = iblock;
  int idist = itot;
  int jdist = 1;
  fftw_r2r_kind kindf[] = {FFTW_R2HC};
  fftw_r2r_kind kindb[] = {FFTW_HC2R};
  iplanf = fftw_plan_many_r2r(1, ni, jmax, fftini, ni, istride, idist,
                              fftouti, ni, istride, idist, kindf, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_many_r2r(1, ni, jmax, fftini, ni, istride, idist,
                              fftouti, ni, istride, idist, kindb, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_many_r2r(1, nj, iblock, fftinj, nj, jstride, jdist,
                              fftoutj, nj, jstride, jdist, kindf, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_many_r2r(1, nj, iblock, fftinj, nj, jstride, jdist,
                              fftoutj, nj, jstride, jdist, kindb, FFTW_EXHAUSTIVE);

  fftwplan = true;

  fftw_forget_wisdom();

  return 0;
}

int cgrid::savefield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename, double offset)
{
  FILE *pFile;
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
    return 1;

  int ijk,jj,kk;

  jj = icells;
  kk = icells*jcells;

  // first, add the offset to the data
  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      for(int i=istart; i<iend; i++)
      {
        ijk = i + j*jj + k*kk;
        tmp1[ijk] = data[ijk] + offset;
      }

  // second, save the data to disk
  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fwrite(&tmp1[ijk], sizeof(double), imax, pFile);
      }

  fclose(pFile);

  return 0;
}

int cgrid::loadfield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename, double offset)
{
  FILE *pFile;
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
    return 1;

  int ijk,jj,kk;

  jj = icells;
  kk = icells*jcells;

  // first, load the data from disk
  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fread(&tmp1[ijk], sizeof(double), imax, pFile);
      }

  fclose(pFile);

  // second, remove the offset
  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      for(int i=istart; i<iend; i++)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] = tmp1[ijk] - offset;
      }

  return 0;
}

int cgrid::fftforward(double * restrict data,   double * restrict tmp1,
                      double * restrict fftini, double * restrict fftouti,
                      double * restrict fftinj, double * restrict fftoutj)
{
  int ij,ijk,jj,kk;

  jj = itot;
  kk = itot*jmax;

  // process the fourier transforms slice by slice
  for(int k=0; k<kblock; k++)
  {
#pragma ivdep
    for(int n=0; n<itot*jmax; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      fftini[ij] = data[ijk];
    }

    fftw_execute(iplanf);

#pragma ivdep
    for(int n=0; n<itot*jmax; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      data[ijk] = fftouti[ij];
    }
  }

  jj = iblock;
  kk = iblock*jtot;

  // do the second fourier transform
  for(int k=0; k<kblock; k++)
  {
#pragma ivdep
    for(int n=0; n<iblock*jtot; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      fftinj[ij] = data[ijk];
    }

    fftw_execute(jplanf);

#pragma ivdep
    for(int n=0; n<iblock*jtot; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      // shift to use p in pressure solver
      data[ijk] = fftoutj[ij];
    }
  }

  return 0;
}

int cgrid::fftbackward(double * restrict data,   double * restrict tmp1,
                       double * restrict fftini, double * restrict fftouti,
                       double * restrict fftinj, double * restrict fftoutj)
{
  int ij,ijk,jj,kk;

  jj = iblock;
  kk = iblock*jtot;

  // transform the second transform back
  for(int k=0; k<kblock; k++)
  {
#pragma ivdep
    for(int n=0; n<iblock*jtot; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      fftinj[ij] = data[ijk];
    }

    fftw_execute(jplanb);

#pragma ivdep
    for(int n=0; n<iblock*jtot; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      data[ijk] = fftoutj[ij] / jtot;
    }
  }

  jj = itot;
  kk = itot*jmax;

  // transform the first transform back
  for(int k=0; k<kblock; k++)
  {
#pragma ivdep
    for(int n=0; n<itot*jmax; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      fftini[ij] = data[ijk];
    }

    fftw_execute(iplanb);

#pragma ivdep
    for(int n=0; n<itot*jmax; n++)
    {
      ij  = n;
      ijk = n + k*kk;
      // swap array here to avoid unnecessary 3d loop
      tmp1[ijk] = fftouti[ij] / itot;
    }
  }

  return 0;
}
int cgrid::savexzslice(double * restrict data, double * restrict tmp, char *filename, int jslice)
{
  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax;

  int count = imax*kmax;

  for(int k=0; k<kmax; k++)
#pragma ivdep
    for(int i=0; i<imax; i++)
    {
      // take the modulus of jslice and jmax to have the right offset within proc
      ijk  = i+igc + (jslice+jgc)*jj + (k+kgc)*kk;
      ijkb = i + k*kkb;
      tmp[ijkb] = data[ijk];
    }

  FILE *pFile;
  pFile = fopen(filename, "wb");
  if(pFile == NULL)
    return 1;

  fwrite(tmp, sizeof(double), count, pFile);
  fclose(pFile);

  return 0;
}

int cgrid::savexyslice(double * restrict data, double * restrict tmp, char *filename, int kslice)
{
  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;

  int count = imax*jmax;

  for(int j=0; j<jmax; j++)
#pragma ivdep
    for(int i=0; i<imax; i++)
    {
      // take the modulus of jslice and jmax to have the right offset within proc
      ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
      ijkb = i + j*jjb;
      tmp[ijkb] = data[ijk];
    }

  FILE *pFile;
  pFile = fopen(filename, "wb");
  if(pFile == NULL)
    return 1;

  fwrite(tmp, sizeof(double), count, pFile);
  fclose(pFile);

  return 0;
}

int cgrid::loadxyslice(double * restrict data, double * restrict tmp, char *filename, int kslice)
{
  int count = imax*jmax;

  FILE *pFile;
  pFile = fopen(filename, "rb");
  if(pFile == NULL)
    return 1;

  fread(tmp, sizeof(double), count, pFile);
  fclose(pFile);

  // put the data back into a field with ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;

  for(int j=0; j<jmax; j++)
#pragma ivdep
    for(int i=0; i<imax; i++)
    {
      ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
      ijkb = i + j*jjb;
      data[ijk] = tmp[ijkb];
    }

  return 0;
}
#endif
