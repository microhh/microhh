#ifndef PARALLEL
#include <fftw3.h>
#include <cstdio>
#include "grid.h"
#include "defines.h"

// MPI functions
int cgrid::initmpi()
{
  // create the MPI types for the cyclic boundary conditions
  int datacount, datablock, datastride;

  // east west
  datacount  = jcells*kcells;
  datablock  = igc;
  datastride = icells;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &eastwestedge);
  MPI_Type_commit(&eastwestedge);

  // north south
  datacount  = kcells;
  datablock  = icells*jgc;
  datastride = icells*jcells;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &northsouthedge);
  MPI_Type_commit(&northsouthedge);

  // transposez
  datacount = imax*jmax*kblock;
  MPI_Type_contiguous(datacount, MPI_DOUBLE, &transposez);
  MPI_Type_commit(&transposez);

  // transposez iblock/jblock/kblock
  datacount = iblock*jblock*kblock;
  MPI_Type_contiguous(datacount, MPI_DOUBLE, &transposez2);
  MPI_Type_commit(&transposez2);

  // transposex imax
  datacount  = jmax*kblock;
  datablock  = imax;
  datastride = itot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposex);
  MPI_Type_commit(&transposex);

  // transposex iblock
  datacount  = jmax*kblock;
  datablock  = iblock;
  datastride = itot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposex2);
  MPI_Type_commit(&transposex2);

  // transposey
  datacount  = kblock;
  datablock  = iblock*jmax;
  datastride = iblock*jtot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposey);
  MPI_Type_commit(&transposey);

  // transposey2
  datacount  = kblock;
  datablock  = iblock*jblock;
  datastride = iblock*jtot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposey2);
  MPI_Type_commit(&transposey2);

  // file saving and loading, take C-ordering into account
  int totsizei  = itot;
  int subsizei  = imax;
  int substarti = mpi->mpicoordx*imax;
  MPI_Type_create_subarray(1, &totsizei, &subsizei, &substarti, MPI_ORDER_C, MPI_DOUBLE, &subi);
  MPI_Type_commit(&subi);

  int totsizej  = jtot;
  int subsizej  = jmax;
  int substartj = mpi->mpicoordy*jmax;
  MPI_Type_create_subarray(1, &totsizej, &subsizej, &substartj, MPI_ORDER_C, MPI_DOUBLE, &subj);
  MPI_Type_commit(&subj);

  // int totsize [3] = {kmax, jtot, itot};
  // int subsize [3] = {kmax, jmax, imax};
  // int substart[3] = {0, mpi->mpicoordy*jmax, mpi->mpicoordx*imax};
  int totsize [3] = {kmax  , jtot, itot};
  int subsize [3] = {kblock, jmax, itot};
  int substart[3] = {mpi->mpicoordx*kblock, mpi->mpicoordy*jmax, 0};
  MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &subarray);
  MPI_Type_commit(&subarray);

  // allocate the array for the profiles
  profl = new double[kcells];

  mpitypes = true;

  return 0;
} 

int cgrid::exitmpi()
{
  if(mpitypes)
  {
    MPI_Type_free(&eastwestedge);
    MPI_Type_free(&northsouthedge);
    MPI_Type_free(&transposez);
    MPI_Type_free(&transposez2);
    MPI_Type_free(&transposex);
    MPI_Type_free(&transposex2);
    MPI_Type_free(&transposey);
    MPI_Type_free(&transposey2);
    MPI_Type_free(&subi);
    MPI_Type_free(&subj);
    MPI_Type_free(&subarray);

    delete[] profl;
  }

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
        for(int i=istart; i<iend; i++)
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
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  if(mpi->mpiid == 0) std::printf("Saving \"%s\"\n", filename);

  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%s\" cannot be written\n", filename);
    return 1;
  }

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  if(mpi->mpicoordy == 0)
    MPI_File_write(fh, &x[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  if(mpi->mpicoordy == 0)
    MPI_File_write(fh, &xh[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  if(mpi->mpicoordx == 0)
    MPI_File_write(fh, &y[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += jtot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  if(mpi->mpicoordx == 0)
    MPI_File_write(fh, &yh[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);

  MPI_File_sync(fh);
  if(MPI_File_close(&fh))
    return 1;

  if(mpi->mpiid == 0)
  {
    FILE *pFile;
    pFile = fopen(filename, "ab");
    fwrite(&z [kstart], sizeof(double), kmax, pFile);
    fwrite(&zh[kstart], sizeof(double), kmax, pFile);
    fclose(pFile);
  }

  // SAVE THE FFTW PLAN
  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", "fftwplan", 0);

    std::printf("Saving \"%s\"\n", filename);

    int n = fftw_export_wisdom_to_filename(filename);
    if(n == 0)
    {
      std::printf("ERROR \"%s\" cannot be saved\n", filename);
      return 1;
    }
  }

  return 0;
}

int cgrid::load()
{
  // LOAD THE GRID
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  if(mpi->mpiid == 0) std::printf("Loading \"%s\"\n", filename);

  FILE *pFile;
  if(mpi->mpiid == 0)
  {
    pFile = fopen(filename, "rb");
    int n = (2*itot+2*jtot)*sizeof(double);
    fseek(pFile, n, SEEK_SET);
    fread(&z [kstart], sizeof(double), kmax, pFile);
    fread(&zh[kstart], sizeof(double), kmax, pFile);
    fclose(pFile);
  }

  mpi->broadcast(&z [kstart], kmax);
  mpi->broadcast(&zh[kstart], kmax);

  // calculate the missing coordinates
  calculate();

  // LOAD THE FFTW PLAN
  std::sprintf(filename, "%s.%07d", "fftwplan", 0);

  if(mpi->mpiid == 0)
    std::printf("Loading \"%s\"\n", filename);

  int n = fftw_import_wisdom_from_filename(filename);
  if(n == 0)
  {
    if(mpi->mpiid == 0)
      std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  fftw_forget_wisdom();

  return 0;
}

int cgrid::savefield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename)
{
  FILE *pFile;
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
    return 1;

  int ijk,jj,kk;

  jj = icells;
  kk = icells*jcells;

  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fwrite(&data[ijk], sizeof(double), imax, pFile);
      }

  fclose(pFile);

  return 0;
}

int cgrid::loadfield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename)
{
  FILE *pFile;
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
    return 1;

  int ijk,jj,kk;

  jj = icells;
  kk = icells*jcells;

  for(int k=kstart; k<kend; k++)
    for(int j=jstart; j<jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fread(&data[ijk], sizeof(double), imax, pFile);
      }

  fclose(pFile);

  return 0;
}

int cgrid::fftforward(double * restrict data,   double * restrict tmp1,
                      double * restrict fftini, double * restrict fftouti,
                      double * restrict fftinj, double * restrict fftoutj)
{
  int ijk,jj,kk;

  // transpose the pressure field
  transposezx(tmp1,data);

  jj = itot;
  kk = itot*jmax;

  // do the first fourier transform
  for(int k=0; k<kblock; k++)
    for(int j=0; j<jmax; j++)
    {
#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        fftini[i] = tmp1[ijk];
      }

      fftw_execute(iplanf);

#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        tmp1[ijk] = fftouti[i];
      }
    }

  // transpose again
  transposexy(data,tmp1);

  jj = iblock;
  kk = iblock*jtot;

  // do the second fourier transform
  for(int k=0; k<kblock; k++)
    for(int i=0; i<iblock; i++)
    {
      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        fftinj[j] = data[ijk];
      }

      fftw_execute(jplanf);

      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        // shift to use p in pressure solver
        tmp1[ijk] = fftoutj[j];
      }
    }

  // transpose back to original orientation
  transposeyz(data,tmp1);

  return 0;
}

int cgrid::fftbackward(double * restrict data,   double * restrict tmp1,
                       double * restrict fftini, double * restrict fftouti,
                       double * restrict fftinj, double * restrict fftoutj)
{
  int ijk,jj,kk;

  // transpose back to y
  transposezy(tmp1, data);

  jj = iblock;
  kk = iblock*jtot;

  // transform the second transform back
  for(int k=0; k<kblock; k++)
    for(int i=0; i<iblock; i++)
    {
      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        fftinj[j] = tmp1[ijk];
      }

      fftw_execute(jplanb);

      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        data[ijk] = fftoutj[j] / jtot;
      }
    }

  // transpose back to x
  transposeyx(tmp1, data);

  jj = itot;
  kk = itot*jmax;

  // transform the first transform back
  for(int k=0; k<kblock; k++)
    for(int j=0; j<jmax; j++)
    {
#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        fftini[i] = tmp1[ijk];
      }

      fftw_execute(iplanb);

#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        // swap array here to avoid unncessary 3d loop
        data[ijk] = fftouti[i] / itot;
      }
    }

  // and transpose back...
  transposexz(tmp1, data);

  jj = imax;
  kk = imax*jmax;

  int ijkp,jjp,kkp1;
  jjp  = icells;
  kkp1 = icells*jcells;

  // put the pressure back onto the original grid including ghost cells
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp1;
        ijk  = i + j*jj + k*kk;
        data[ijkp] = tmp1[ijk];
      }

  return 0;
}
#endif
