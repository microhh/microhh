#ifdef PARALLEL
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

  // east west 2d
  datacount  = jcells;
  datablock  = igc;
  datastride = icells;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &eastwestedge2d);
  MPI_Type_commit(&eastwestedge2d);

  // north south 2d
  datacount  = 1;
  datablock  = icells*jgc;
  datastride = icells*jcells;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &northsouthedge2d);
  MPI_Type_commit(&northsouthedge2d);

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

  // the lines below describe the array in case transposes are not used before saving
  // int totsize [3] = {kmax, jtot, itot};
  // int subsize [3] = {kmax, jmax, imax};
  // int substart[3] = {0, mpi->mpicoordy*jmax, mpi->mpicoordx*imax};
  int totsize [3] = {kmax  , jtot, itot};
  int subsize [3] = {kblock, jmax, itot};
  int substart[3] = {mpi->mpicoordx*kblock, mpi->mpicoordy*jmax, 0};
  MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &subarray);
  MPI_Type_commit(&subarray);

  // save mpitype for a xz-slice for cross section processing
  int totxzsize [2] = {kmax, itot};
  int subxzsize [2] = {kmax, imax};
  int subxzstart[2] = {0, mpi->mpicoordx*imax};
  MPI_Type_create_subarray(2, totxzsize, subxzsize, subxzstart, MPI_ORDER_C, MPI_DOUBLE, &subxzslice);
  MPI_Type_commit(&subxzslice);

  // save mpitype for a xy-slice for cross section processing
  int totxysize [2] = {jtot, itot};
  int subxysize [2] = {jmax, imax};
  int subxystart[2] = {mpi->mpicoordy*jmax, mpi->mpicoordx*imax};
  MPI_Type_create_subarray(2, totxysize, subxysize, subxystart, MPI_ORDER_C, MPI_DOUBLE, &subxyslice);
  MPI_Type_commit(&subxyslice);

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
    MPI_Type_free(&subxzslice);

    delete[] profl;
  }

  return 0;
}

int cgrid::boundary_cyclic(double * restrict data)
{
  int ncount = 1;

  // communicate east-west edges
  int eastout = iend-igc;
  int westin  = 0;
  int westout = istart;
  int eastin  = iend;

  // communicate north-south edges
  int northout = (jend-jgc)*icells;
  int southin  = 0;
  int southout = jstart*icells;
  int northin  = jend  *icells;

  // first, send and receive the ghost cells in east-west direction
  MPI_Isend(&data[eastout], ncount, eastwestedge, mpi->neast, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Irecv(&data[westin], ncount, eastwestedge, mpi->nwest, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Isend(&data[westout], ncount, eastwestedge, mpi->nwest, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Irecv(&data[eastin], ncount, eastwestedge, mpi->neast, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  // wait here for the mpi to have correct values in the corners of the cells
  mpi->waitall();

  // if the run is 3D, apply the BCs
  if(jtot > 1)
  {
    // second, send and receive the ghost cells in the north-south direction
    MPI_Isend(&data[northout], ncount, northsouthedge, mpi->nnorth, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&data[southin], ncount, northsouthedge, mpi->nsouth, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Isend(&data[southout], ncount, northsouthedge, mpi->nsouth, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&data[northin], ncount, northsouthedge, mpi->nnorth, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    mpi->waitall();
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

int cgrid::boundary_cyclic2d(double * restrict data)
{
  int ncount = 1;

  // communicate east-west edges
  int eastout = iend-igc;
  int westin  = 0;
  int westout = istart;
  int eastin  = iend;

  // communicate north-south edges
  int northout = (jend-jgc)*icells;
  int southin  = 0;
  int southout = jstart*icells;
  int northin  = jend  *icells;

  // first, send and receive the ghost cells in east-west direction
  MPI_Isend(&data[eastout], ncount, eastwestedge2d, mpi->neast, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Irecv(&data[westin], ncount, eastwestedge2d, mpi->nwest, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Isend(&data[westout], ncount, eastwestedge2d, mpi->nwest, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Irecv(&data[eastin], ncount, eastwestedge2d, mpi->neast, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  // wait here for the mpi to have correct values in the corners of the cells
  mpi->waitall();

  // if the run is 3D, apply the BCs
  if(jtot > 1)
  {
    // second, send and receive the ghost cells in the north-south direction
    MPI_Isend(&data[northout], ncount, northsouthedge2d, mpi->nnorth, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&data[southin], ncount, northsouthedge2d, mpi->nsouth, 1, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Isend(&data[southout], ncount, northsouthedge2d, mpi->nsouth, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&data[northin], ncount, northsouthedge2d, mpi->nnorth, 2, mpi->commxy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    mpi->waitall();
  }
  // in case of 2D, fill all the ghost cells with the current value
  else
  {
    // 2d essential variables
    int ijref,ijnorth,ijsouth,jj;

    jj = icells;

    for(int j=0; j<jgc; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
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
  int ijks, ijkr;
  int ncount = 1;
  int tag = 1;

  int jj = imax;
  int kk = imax*jmax;

  for(int n=0; n<mpi->npx; n++)
  {
    // determine where to fetch the data and where to store it
    ijks = n*kblock*kk;
    ijkr = n*jj;

    // send and receive the data
    MPI_Isend(&as[ijks], ncount, transposez, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&ar[ijkr], ncount, transposex, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }
  mpi->waitall();

  return 0;
}

int cgrid::transposexz(double * restrict ar, double * restrict as)
{
  int ijks, ijkr;
  int ncount = 1;
  int tag = 1;

  int jj = imax;
  int kk = imax*jmax;

  for(int n=0; n<mpi->npx; n++)
  {
    // determine where to fetch the data and where to store it
    ijks = n*jj;
    ijkr = n*kblock*kk;

    // send and receive the data
    MPI_Isend(&as[ijks], ncount, transposex, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&ar[ijkr], ncount, transposez, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }
  mpi->waitall();

  return 0;
}

int cgrid::transposexy(double * restrict ar, double * restrict as)
{
  int ijks, ijkr;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jmax;

  for(int n=0; n<mpi->npy; n++)
  {
    // determine where to fetch the data and where to store it
    ijks = n*jj;
    ijkr = n*kk;

    // send and receive the data
    MPI_Isend(&as[ijks], ncount, transposex2, n, tag, mpi->commy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&ar[ijkr], ncount, transposey , n, tag, mpi->commy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }
  mpi->waitall();

  return 0;
}

int cgrid::transposeyx(double * restrict ar, double * restrict as)
{
  int ijks, ijkr;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jmax;

  for(int n=0; n<mpi->npy; n++)
  {
    // determine where to fetch the data and where to store it
    ijks = n*kk;
    ijkr = n*jj;

    // send and receive the data
    MPI_Isend(&as[ijks], ncount, transposey , n, tag, mpi->commy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&ar[ijkr], ncount, transposex2, n, tag, mpi->commy, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }
  mpi->waitall();

  return 0;
}

int cgrid::transposeyz(double * restrict ar, double * restrict as)
{
  int ijks,ijkr;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jblock;

  for(int n=0; n<mpi->npx; n++)
  {
    // determine where to fetch the data and where to store it
    ijks = n*jblock*jj;
    ijkr = n*kblock*kk;

    // send and receive the data
    MPI_Isend(&as[ijks], ncount, transposey2, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&ar[ijkr], ncount, transposez2, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }
  mpi->waitall();
 
  return 0;
}

int cgrid::transposezy(double * restrict ar, double * restrict as)
{
  int ijks,ijkr;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jblock;

  for(int n=0; n<mpi->npx; n++)
  {
    // determine where to fetch the data and where to store it
    ijks = n*kblock*kk;
    ijkr = n*jblock*jj;

    // send and receive the data
    MPI_Isend(&as[ijks], ncount, transposez2, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&ar[ijkr], ncount, transposey2, n, tag, mpi->commx, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }
  mpi->waitall();

  return 0;
}

int cgrid::getmax(double *var)
{
  double varl = *var;
  MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_MAX, mpi->commxy);

  return 0;
}

int cgrid::getsum(double *var)
{
  double varl = *var;
  MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_SUM, mpi->commxy);

  return 0;
}

int cgrid::getprof(double *prof, int kcellsin)
{
  for(int k=0; k<kcellsin; k++)
    profl[k] = prof[k] / mpi->nprocs;

  MPI_Allreduce(profl, prof, kcellsin, MPI_DOUBLE, MPI_SUM, mpi->commxy);

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

  fftwplan = true;

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
    if(pFile == NULL)
    {
      std::printf("ERROR \"%s\" does not exist\n", filename);
      return 1;
    }
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

  fftwplan = true;

  fftw_forget_wisdom();

  return 0;
}

int cgrid::savefield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename)
{
  // save the data in transposed order to have large chunks of contiguous disk space
  // MPI-IO is not stable on Juqueen and supermuc otherwise

  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax*jmax;

  int count = imax*jmax*kmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        ijkb = i + j*jjb + k*kkb;
        tmp1[ijkb] = data[ijk];
      }

  transposezx(tmp2, tmp1);

  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  if(MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL))
    return 1;

  if(MPI_File_write_all(fh, tmp2, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
    return 1;

  if(MPI_File_close(&fh))
    return 1;

  return 0;
}

int cgrid::loadfield3d(double * restrict data, double * restrict tmp1, double * restrict tmp2, char *filename)
{
  // save the data in transposed order to have large chunks of contiguous disk space
  // MPI-IO is not stable on Juqueen and supermuc otherwise

  // read the file
  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";
  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL);

  // extract the data from the 3d field without the ghost cells
  int count = imax*jmax*kmax;

  if(MPI_File_read_all(fh, tmp1, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
    return 1;

  if(MPI_File_close(&fh))
    return 1;

  // transpose the data back
  transposexz(tmp2, tmp1);

  int ijk,jj,kk;
  int ijkb,jjb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax*jmax;

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        ijkb = i + j*jjb + k*kkb;
        data[ijk] = tmp2[ijkb];
      }

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

  return 0;
}

int cgrid::savexzslice(double * restrict data, double * restrict tmp, int jslice, char *filename)
{
  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,kkb;

  jj  = icells;
  kk  = icells*jcells;
  kkb = imax;

  int count = imax*kmax;

  for(int k=0; k<kmax; k++)
#pragma ivdep
    for(int i=0; i<imax; i++)
    {
      // take the modulus of jslice and jmax to have the right offset within proc
      ijk  = i+igc + ((jslice%jmax)+jgc)*jj + (k+kgc)*kk;
      ijkb = i + k*kkb;
      tmp[ijkb] = data[ijk];
    }

  if(mpi->mpicoordy == jslice/jmax)
  {
    MPI_File fh;
    if(MPI_File_open(mpi->commx, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
      return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if(MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subxzslice, name, MPI_INFO_NULL))
      return 1;

    // only write at the procs that contain the slice
    if(MPI_File_write_all(fh, tmp, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
      return 1;

    MPI_File_sync(fh);

    if(MPI_File_close(&fh))
      return 1;
  }

  MPI_Barrier(mpi->commxy);

  return 0;
}

int cgrid::savexyslice(double * restrict data, double * restrict tmp, int kslice, char *filename)
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

  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  if(MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subxyslice, name, MPI_INFO_NULL))
    return 1;

  // only write at the procs that contain the slice
  if(MPI_File_write_all(fh, tmp, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
    return 1;

  MPI_File_sync(fh);

  if(MPI_File_close(&fh))
    return 1;

  MPI_Barrier(mpi->commxy);

  return 0;
}
#endif
