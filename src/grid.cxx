#include <cstdio>
#include <cmath>
#include "grid.h"
#include "input.h"
#include "defines.h"

// build the grid
cgrid::cgrid(cmpi *mpiin)
{
  // std::printf("Creating instance of object grid\n");
  mpi = mpiin;

  allocated = false;
}

cgrid::~cgrid()
{
  if(allocated)
  { 
    delete[] x;
    delete[] xh;
    delete[] y;
    delete[] yh;
    delete[] z;
    delete[] zh;
    delete[] dz;
    delete[] dzh;
    delete[] dzi;
    delete[] dzhi;
    delete[] dzi4;
    delete[] dzhi4;
  }

  // std::printf("Destroying instance of object grid\n");
}

int cgrid::readinifile(cinput *inputin)
{
  /*// setup Taylor-Green vortex
  xsize = 1.;
  ysize = 1.;
  zsize = 0.5;

  itot  = 64;
  jtot  = 8;
  ktot  = 32;
  // end setup Taylor-Green vortex*/

  int n = 0;

  n += inputin->getItem(&xsize, "grid", "xsize");
  n += inputin->getItem(&ysize, "grid", "ysize");
  n += inputin->getItem(&zsize, "grid", "zsize");

  n += inputin->getItem(&itot, "grid", "itot");
  n += inputin->getItem(&jtot, "grid", "jtot");
  n += inputin->getItem(&ktot, "grid", "ktot");

  if(n > 0)
    return 1;
  
  igc = 3;
  jgc = 3;
  kgc = 3;

  return 0;
}

int cgrid::init()
{
  // check whether the grid fits the processor configuration
  if(itot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npx = %d\n", itot, mpi->npx);
    return 1;
  }
  if(itot % mpi->npy != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npy = %d\n", itot, mpi->npy);
    return 1;
  }
  // check this one only when npy > 1, since the transpose in that direction only happens then
  //if(jtot % mpi->npx != 0 && mpi->npy > 1)
  if(jtot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npx = %d\n", jtot, mpi->npx);
    return 1;
  }
  if(jtot % mpi->npy != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npy = %d\n", jtot, mpi->npy);
    return 1;
  }
  if(ktot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR ktot = %d is not a multiple of npx = %d\n", ktot, mpi->npx);
    return 1;
  }

  imax   = itot / mpi->npx;
  jmax   = jtot / mpi->npy;
  kmax   = ktot;

  iblock = itot / mpi->npy;
  jblock = jtot / mpi->npx;
  kblock = ktot / mpi->npx;

  icells = (imax+2*igc);
  jcells = (jmax+2*jgc);
  kcells = (kmax+2*kgc);
  ncells = (imax+2*igc)*(jmax+2*jgc)*(kmax+2*kgc);

  istart = igc;
  jstart = jgc;
  kstart = kgc;

  iend   = imax + igc;
  jend   = jmax + jgc;
  kend   = kmax + kgc;

  x     = new double[imax+2*igc];
  xh    = new double[imax+2*igc];
  y     = new double[jmax+2*jgc];
  yh    = new double[jmax+2*jgc];
  z     = new double[kmax+2*kgc];
  zh    = new double[kmax+2*kgc];
  dz    = new double[kmax+2*kgc];
  dzh   = new double[kmax+2*kgc];
  dzi   = new double[kmax+2*kgc];
  dzhi  = new double[kmax+2*kgc];
  dzi4  = new double[kmax+2*kgc];
  dzhi4 = new double[kmax+2*kgc];

  allocated = true;

  // initialize the communication functions
  initmpi();

  return 0;
}

int cgrid::create(cinput *inputin)
{
  if(inputin->getProf(&z[kstart], "z", kmax))
    return 1;

  calculate();

  return 0;
}

int cgrid::calculate()
{
  int i,j,k;

  dx = xsize / itot;
  dy = ysize / jtot;

  double xoff = mpi->mpicoordx * xsize / mpi->npx;
  double yoff = mpi->mpicoordy * ysize / mpi->npy;

  // calculate the x and y coordinates
  for(i=0; i<icells; i++)
  {
    x [i] = 0.5*dx + (i-igc)*dx + xoff;
    xh[i] = (i-igc)*dx + xoff;
  }

  for(j=0; j<jcells; j++)
  {
    y [j] = 0.5*dy + (j-jgc)*dy + yoff;
    yh[j] = (j-jgc)*dy + yoff;
  }

  // calculate the height of the ghost cell
  z[kstart-1] = -z[kstart  ];
  z[kstart-2] = -z[kstart+1];
  z[kstart-3] = -z[kstart+2];

  z[kend  ] = 2.*zsize - z[kend-1];
  z[kend+1] = 2.*zsize - z[kend-2];
  z[kend+2] = 2.*zsize - z[kend-3];

  // calculate the half levels according to the numerical scheme
  // zh[0] = -999.;
  // zh[kstart-1] = interp4biasbot(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1]);
  zh[kstart  ] = 0.;
  for(k=kstart+1; k<kend; k++)
    zh[k] = interp4(z[k-2], z[k-1], z[k], z[k+1]);
  zh[kend] = zsize;
  // zh[kend+1] = interp4biastop(z[kend-2], z[kend-1], z[kend], z[kend+1]);
  //
  zh[kstart-1] = -zh[kstart+1];
  zh[kstart-2] = -zh[kstart+2];
  zh[kstart-3] = -zh[kstart+3];

  zh[kend+1] = 2.*zsize - zh[kend-1];
  zh[kend+2] = 2.*zsize - zh[kend-2];

  // compute the height of the grid cells
  for(k=1; k<kcells; k++)
  {
    dzh [k] = z[k] - z[k-1];
    dzhi[k] = 1./dzh[k];
  }
  dzh [kstart-3] = dzh [kstart+3];
  dzhi[kstart-3] = dzhi[kstart+3];

  // compute the height of the grid cells
  for(k=1; k<kcells-1; k++)
  {
    dz [k] = zh[k+1] - zh[k];
    dzi[k] = 1./dz[k];
  }
  dz [kstart-3] = dz [kstart+2];
  dzi[kstart-3] = dzi[kstart+2];
  dz [kend+2] = dz [kend-2];
  dzi[kend+2] = dzi[kend-2];

  /*
  // calculate the inverse gradients for the 4th order scheme
  dzi4 [kstart] = 1./grad4xbiasbot(zh[kstart  ], zh[kstart+1], zh[kstart+2], zh[kstart+3]);
  dzhi4[kstart] = 1./grad4xbiasbot(z [kstart-1], z [kstart  ], z [kstart+1], z [kstart+2]);
  for(k=kstart+1; k<kend-1; k++)
  {
    dzi4 [k] = 1./grad4x(zh[k-1], zh[k  ], zh[k+1], zh[k+2]);
    dzhi4[k] = 1./grad4x(z [k-2], z [k-1], z [k  ], z [k+1]);
  }
  dzi4 [kend-1] = 1./grad4xbiastop(zh[kend-3], zh[kend-2], zh[kend-1], zh[kend]);

  dzhi4[kend-1] = 1./grad4x       (z[kend-3], z[kend-2], z[kend-1], z[kend]);
  dzhi4[kend  ] = 1./grad4xbiastop(z[kend-3], z[kend-2], z[kend-1], z[kend]);
  */
  
  // calculate the inverse gradients for the 4th order scheme
  // dzi4 [0] = -999.;
  // dzhi4[0] = -999.;
  // dzi4 [kstart-1] = 1./grad4xbiasbot(zh[kstart-1], zh[kstart  ], zh[kstart+1], zh[kstart+2]);
  // dzhi4[kstart-1] = 1./grad4xbiasbot(z [kstart-2], z [kstart-1], z [kstart  ], z [kstart+1]);
  //
  for(k=kstart; k<kend; k++)
  {
    dzi4 [k] = 1./grad4x(zh[k-1], zh[k  ], zh[k+1], zh[k+2]);
    dzhi4[k] = 1./grad4x(z [k-2], z [k-1], z [k  ], z [k+1]);
  }
  dzhi4[kend  ] = 1./grad4x(z [kend-2], z [kend-1], z [kend], z [kend+1]);

  // bc's
  dzi4 [kstart-3] = dzi4 [kstart+2];
  dzhi4[kstart-3] = dzhi4[kstart+3];
  dzi4 [kstart-2] = dzi4 [kstart+1];
  dzhi4[kstart-2] = dzhi4[kstart+2];
  dzi4 [kstart-1] = dzi4 [kstart  ];
  dzhi4[kstart-1] = dzhi4[kstart+1];

  dzi4 [kend  ] = dzi4 [kend-1];
  dzi4 [kend+1] = dzi4 [kend-2];
  dzhi4[kend+1] = dzhi4[kend-1];
  dzi4 [kend+2] = dzi4 [kend-3];
  dzhi4[kend+2] = dzhi4[kend-2];

  return 0;
}

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

  int totsize [3] = {kmax, jtot, itot};
  int subsize [3] = {kmax, jmax, imax};
  int substart[3] = {0, mpi->mpicoordy*jmax, mpi->mpicoordx*imax};
  MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &subarray);
  MPI_Type_commit(&subarray);

  return 0;
} 

#ifdef MPISUBCOMM
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

  mpi->waitall();

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

// IO functions
int cgrid::save()
{
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  if(mpi->mpiid == 0) std::printf("Saving \"%s\"\n", filename);

  /*FILE *pFile;
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" cannot be written\n", filename);
    return 1;
  }
  else
    std::printf("Saving \"%s\"\n", filename);

  fwrite(&x [istart], sizeof(double), imax, pFile);
  fwrite(&xh[istart], sizeof(double), imax, pFile);
  fwrite(&y [jstart], sizeof(double), jmax, pFile);
  fwrite(&yh[jstart], sizeof(double), jmax, pFile);
  fwrite(&z [kstart], sizeof(double), kmax, pFile);
  fwrite(&zh[kstart], sizeof(double), kmax, pFile);
  fclose(pFile);*/

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

  return 0;
}

int cgrid::load()
{
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  if(mpi->mpiid == 0) std::printf("Loading \"%s\"\n", filename);

  /*
  FILE *pFile;
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }
  else
    std::printf("Loading \"%s\"\n", filename);

  fread(&x [istart], sizeof(double), imax, pFile);
  fread(&xh[istart], sizeof(double), imax, pFile);
  fread(&y [jstart], sizeof(double), jmax, pFile);
  fread(&yh[jstart], sizeof(double), jmax, pFile);
  fread(&z [kstart], sizeof(double), kmax, pFile);
  fread(&zh[kstart], sizeof(double), kmax, pFile);
  fclose(pFile);
  */

  // DISABLE READING OF THE HORIZONTAL DATA, IT SLOWS DOWN THE INIT UNNECESSARILY
  /*
  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%s\" cannot be loaded\n", filename);
    return 1;
  }

  // select noncontiguous part of 3d array to store the selected data

  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  // if(mpi->mpiid / mpi->npx == 0)
  MPI_File_read_all(fh, &x[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  // if(mpi->mpiid / mpi->npx == 0)
  MPI_File_read_all(fh, &xh[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  // if(mpi->mpiid % mpi->npx == 0)
  MPI_File_read_all(fh, &y[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += jtot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  // if(mpi->mpiid % mpi->npx == 0)
  MPI_File_read_all(fh, &yh[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);

  MPI_File_sync(fh);
  if(MPI_File_close(&fh))
    return 1;
   */

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
  MPI_Bcast(&z [kstart], kmax, MPI_DOUBLE, 0, mpi->commxy);
  MPI_Bcast(&zh[kstart], kmax, MPI_DOUBLE, 0, mpi->commxy);

  // calculate the missing coordinates
  calculate();

  return 0;
}

int cgrid::savefield3d(double * restrict data, char *filename)
{
  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  if(MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL))
    return 1;

  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;
  // int igc,jgc,kgc;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax*jmax;
  // igc = igc;
  // jgc = jgc;
  // kgc = kgc;

  int count = imax*jmax*kmax;

  double *buffer = new double[count];

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        ijkb = i + j*jjb + k*kkb;
        buffer[ijkb] = data[ijk];
      }

  if(MPI_File_write_all(fh, buffer, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
  {
    delete[] buffer;
    return 1;
  }

  if(MPI_File_close(&fh))
  {
    delete[] buffer;
    return 1;
  }

  delete[] buffer;

  return 0;
}

int cgrid::loadfield3d(double *data, char *filename)
{
  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";
  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL);

  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;
  // int igc,jgc,kgc;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax*jmax;
  // igc = igc;
  // jgc = jgc;
  // kgc = kgc;

  int count = imax*jmax*kmax;
  double *buffer = new double[count];

  if(MPI_File_read_all(fh, buffer, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
  {
    delete[] buffer;
    return 1;
  }

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        ijkb = i + j*jjb + k*kkb;
        data[ijk] = buffer[ijkb];
      }

  if(MPI_File_close(&fh))
  {
    delete[] buffer;
    return 1;
  }

  delete[] buffer;

  return 0;
}

#else
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
  MPI_Isend(&data[eastout], ncount, eastwestedge, mpi->neast, 1, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Irecv(&data[westin], ncount, eastwestedge, mpi->nwest, 1, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Isend(&data[westout], ncount, eastwestedge, mpi->nwest, 2, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;
  MPI_Irecv(&data[eastin], ncount, eastwestedge, mpi->neast, 2, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
  mpi->reqsn++;

  // if the run is 3D, apply the BCs
  if(jtot > 1)
  {
    // second, send and receive the ghost cells in the north-south direction
    MPI_Isend(&data[northout], ncount, northsouthedge, mpi->nnorth, 1, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&data[southin], ncount, northsouthedge, mpi->nsouth, 1, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Isend(&data[southout], ncount, northsouthedge, mpi->nsouth, 2, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
    MPI_Irecv(&data[northin], ncount, northsouthedge, mpi->nnorth, 2, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
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

  mpi->waitall();

  return 0;
}

int cgrid::transposezx(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;
  int tag = 1;

  int jj = imax;
  int kk = imax*jmax;

  for(int k=0; k<mpi->npx; k++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + k;

    // determine where to fetch the send data
    int ijks = k*kblock*kk;

    // send the block
    MPI_Isend(&as[ijks], ncount, transposez, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;

    // determine where to store the receive data
    int ijkr = k*jj;

    // receive the block
    MPI_Irecv(&ar[ijkr], ncount, transposex, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }

  mpi->waitall();

  return 0;
}

int cgrid::transposexz(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;
  int tag = 1;

  int jj = imax;
  int kk = imax*jmax;

  for(int i=0; i<mpi->npx; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + i;

    // determine where to fetch the send data
    int ijks = i*jj;

    // send the block
    MPI_Isend(&as[ijks], ncount, transposex, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;

    // determine where to store the receive data
    int ijkr = i*kblock*kk;

    // receive the block
    MPI_Irecv(&ar[ijkr], ncount, transposez, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }

  mpi->waitall();

  return 0;
}

int cgrid::transposexy(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jmax;

  for(int i=0; i<mpi->npy; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid % mpi->npx + i * mpi->npx;

    // determine where to fetch the send data
    int ijks = i*jj;

    // send the block
    MPI_Isend(&as[ijks], ncount, transposex2, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;

    // determine where to store the receive data
    int ijkr = i*kk;

    // receive the block
    MPI_Irecv(&ar[ijkr], ncount, transposey, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }

  mpi->waitall();

  return 0;
}

int cgrid::transposeyx(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jmax;

  for(int i=0; i<mpi->npy; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid % mpi->npx + i * mpi->npx;

    // determine where to fetch the send data
    int ijks = i*kk;

    // send the block
    MPI_Isend(&as[ijks], ncount, transposey, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;

    // determine where to store the receive data
    int ijkr = i*jj;

    // receive the block
    MPI_Irecv(&ar[ijkr], ncount, transposex2, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }

  mpi->waitall();

  return 0;
}

int cgrid::transposeyz(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jblock;

  for(int i=0; i<mpi->npx; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + i;

    // determine where to fetch the send data
    int ijks = i*jblock*jj;

    // send the block, tag it with the height (in kblocks) where it should come
    MPI_Isend(&as[ijks], ncount, transposey2, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;

    // determine where to store the receive data
    int ijkr = i*kblock*kk;

    // and determine what has to be delivered at height i (in kblocks)
    MPI_Irecv(&ar[ijkr], ncount, transposez2, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }

  mpi->waitall();
 
  return 0;
}

int cgrid::transposezy(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;
  int tag = 1;

  int jj = iblock;
  int kk = iblock*jblock;

  for(int k=0; k<mpi->npx; k++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + k;

    // determine where to fetch the send data
    int ijks = k*kblock*kk;

    // send the block
    MPI_Isend(&as[ijks], ncount, transposez2, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;

    // determine where to store the receive data
    int ijkr = k*jblock*jj;

    // receive the block
    MPI_Irecv(&ar[ijkr], ncount, transposey2, nblock, tag, MPI_COMM_WORLD, &mpi->reqs[mpi->reqsn]);
    mpi->reqsn++;
  }

  mpi->waitall();

  return 0;
}

int cgrid::getmax(double *var)
{
  double varl = *var;
  MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return 0;
}

int cgrid::getsum(double *var)
{
  double varl = *var;
  MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return 0;
}

// IO functions
int cgrid::save()
{
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  if(mpi->mpiid == 0) std::printf("Saving \"%s\"\n", filename);

  /*FILE *pFile;
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" cannot be written\n", filename);
    return 1;
  }
  else
    std::printf("Saving \"%s\"\n", filename);

  fwrite(&x [istart], sizeof(double), imax, pFile);
  fwrite(&xh[istart], sizeof(double), imax, pFile);
  fwrite(&y [jstart], sizeof(double), jmax, pFile);
  fwrite(&yh[jstart], sizeof(double), jmax, pFile);
  fwrite(&z [kstart], sizeof(double), kmax, pFile);
  fwrite(&zh[kstart], sizeof(double), kmax, pFile);
  fclose(pFile);*/

  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
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
  MPI_Barrier(MPI_COMM_WORLD);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  if(mpi->mpicoordy == 0)
    MPI_File_write(fh, &xh[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  if(mpi->mpicoordx == 0)
    MPI_File_write(fh, &y[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
  fileoff += jtot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  if(mpi->mpicoordx == 0)
    MPI_File_write(fh, &yh[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);

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

  return 0;
}

int cgrid::load()
{
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  if(mpi->mpiid == 0) std::printf("Loading \"%s\"\n", filename);

  /*
  FILE *pFile;
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }
  else
    std::printf("Loading \"%s\"\n", filename);

  fread(&x [istart], sizeof(double), imax, pFile);
  fread(&xh[istart], sizeof(double), imax, pFile);
  fread(&y [jstart], sizeof(double), jmax, pFile);
  fread(&yh[jstart], sizeof(double), jmax, pFile);
  fread(&z [kstart], sizeof(double), kmax, pFile);
  fread(&zh[kstart], sizeof(double), kmax, pFile);
  fclose(pFile);
  */

  // DISABLE READING OF THE HORIZONTAL DATA, IT SLOWS DOWN THE INIT UNNECESSARILY
  /*
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%s\" cannot be loaded\n", filename);
    return 1;
  }

  // select noncontiguous part of 3d array to store the selected data

  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  // if(mpi->mpiid / mpi->npx == 0)
  MPI_File_read_all(fh, &x[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  // if(mpi->mpiid / mpi->npx == 0)
  MPI_File_read_all(fh, &xh[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  // if(mpi->mpiid % mpi->npx == 0)
  MPI_File_read_all(fh, &y[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
  fileoff += jtot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  // if(mpi->mpiid % mpi->npx == 0)
  MPI_File_read_all(fh, &yh[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_sync(fh);
  if(MPI_File_close(&fh))
    return 1;
   */

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
  MPI_Bcast(&z [kstart], kmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&zh[kstart], kmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // calculate the missing coordinates
  calculate();

  return 0;
}

int cgrid::savefield3d(double * restrict data, char *filename)
{
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  if(MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL))
    return 1;

  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;
  // int igc,jgc,kgc;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax*jmax;
  // igc = igc;
  // jgc = jgc;
  // kgc = kgc;

  int count = imax*jmax*kmax;

  double *buffer = new double[count];

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        ijkb = i + j*jjb + k*kkb;
        buffer[ijkb] = data[ijk];
      }

  if(MPI_File_write_all(fh, buffer, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
  {
    delete[] buffer;
    return 1;
  }

  if(MPI_File_close(&fh))
  {
    delete[] buffer;
    return 1;
  }

  delete[] buffer;

  return 0;
}

int cgrid::loadfield3d(double *data, char *filename)
{
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";
  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL);

  // extract the data from the 3d field without the ghost cells
  int ijk,jj,kk;
  int ijkb,jjb,kkb;
  // int igc,jgc,kgc;

  jj  = icells;
  kk  = icells*jcells;
  jjb = imax;
  kkb = imax*jmax;
  // igc = igc;
  // jgc = jgc;
  // kgc = kgc;

  int count = imax*jmax*kmax;
  double *buffer = new double[count];

  if(MPI_File_read_all(fh, buffer, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
  {
    delete[] buffer;
    return 1;
  }

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        ijkb = i + j*jjb + k*kkb;
        data[ijk] = buffer[ijkb];
      }

  if(MPI_File_close(&fh))
  {
    delete[] buffer;
    return 1;
  }

  delete[] buffer;

  return 0;
}
#endif

inline double cgrid::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

inline double cgrid::interp4biasbot(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*a + (15./16.)*b - (5./16.)*c + (1./16)*d);
}

inline double cgrid::interp4biastop(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*d + (15./16.)*c - (5./16.)*b + (1./16)*a);
}

inline double cgrid::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b)); 
}

inline double cgrid::grad4xbiasbot(const double a, const double b, const double c, const double d)
{
  return (-23.*a + 21.*b + 3.*c - d);
}

inline double cgrid::grad4xbiastop(const double a, const double b, const double c, const double d)
{
  return ( 23.*d - 21.*c - 3.*b + a);
}
