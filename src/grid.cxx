#include <cstdio>
#include <cmath>
#include "grid.h"
#include "input.h"
#include "defines.h"

// build the grid
cgrid::cgrid(cmpi *mpiin)
{
  std::printf("Creating instance of object grid\n");

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
  }

  std::printf("Destroying instance of object grid\n");
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
  
  igc = 1;
  jgc = 1;
  kgc = 1;

  return 0;
}

int cgrid::init()
{
  if(itot % mpi->npx != 0 || itot % mpi->npy != 0)
  {
    std::printf("ERROR itot = %d is not a multiple of npx = %d or npy = %d\n", itot, mpi->npx, mpi->npy);
    return 1;
  }
  if(jtot % mpi->npx != 0 || jtot % mpi->npy != 0)
  {
    std::printf("ERROR jtot = %d is not a multiple of npx = %d or npy = %d\n", jtot, mpi->npx, mpi->npy);
    return 1;
  }
  if(ktot % mpi->npx != 0)
  {
    std::printf("ERROR ktot = %d is not a multiple of npx = %d\n", ktot, mpi->npx);
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

  x    = new double[imax+2*igc];
  xh   = new double[imax+2*igc];
  y    = new double[jmax+2*jgc];
  yh   = new double[jmax+2*jgc];
  z    = new double[kmax+2*kgc];
  zh   = new double[kmax+2*kgc];
  dz   = new double[kmax+2*kgc];
  dzh  = new double[kmax+2*kgc];
  dzi  = new double[kmax+2*kgc];
  dzhi = new double[kmax+2*kgc];

  allocated = true;

  // initialize the communication functions
  initmpi();

  return 0;
}

int cgrid::create()
{
  // create non-equidistant grid
  double alpha = 0.967;
  double eta;
  int k;

  // heights are set according to Moser180 case
  for(k=kstart; k<kend; k++)
  {
    eta  = -1. + 2.*((k-kstart+1) - 0.5) / kmax;
    z[k] = zsize / (2.*alpha) * std::tanh(eta*0.5*(std::log(1.+alpha) - std::log(1.-alpha))) + 0.5*zsize;
  }
  // end Moser180 setup 
  
  /*// uniform height setup
  for(k=kstart; k<kend; k++)
    z[k] = zsize / (2*kmax) + zsize / kmax * (k-kstart);
  // end uniform height setup*/

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

  // calculate the height of the ghost cells
  for(k=0; k<kgc; k++)
  {
    z[kstart-k-1] = -1. * z[kstart+k];
    z[kend  +k  ] = -1. * z[kend-1-k] + 2.*zsize;
  }

  // assume the flux levels are exactly in between the cells
  // compute the flux levels and the distance between them
  for(k=1; k<kcells; k++)
  {
    zh  [k] = 0.5*(z[k] + z[k-1]);
    dzh [k] = z[k] - z[k-1];
    dzhi[k] = 1./dzh[k];
  }

  // set the non-initialized values
  zh  [0] = -zh[2];
  dzh [0] = -999.;
  dzhi[0] = -999.;

  // compute the height of the grid cells
  for(k=kstart; k<kend; k++)
  {
    dz [k] = 0.5*(z[k]-z[k-1]) + 0.5*(z[k+1]-z[k]);
    dzi[k] = 1./dz[k];
  }

  // compute the height of the ghost cells
  for(k=0; k<kgc; k++)
  {
    dz[kstart-k-1]  = dz[kstart+k];
    dz[kend+k]      = dz[kend-k-1];
    dzi[kstart-k-1] = 1./dz[kstart-k-1];
    dzi[kend+k]     = 1./dz[kend+k];
  }

  return 0;
}

int cgrid::save()
{
  char filename[256];
  std::sprintf(filename, "%s.%07d", "grid", 0);
  if(mpi->mpiid == 0)
    std::printf("Saving \"%s\"\n", filename);

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
    if(mpi->mpiid == 0)
      std::printf("ERROR \"%s\" cannot be written\n", filename);
    return 1;
  }

  // select noncontiguous part of 3d array to store the selected data

  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  if(mpi->mpiid / mpi->npx == 0)
    MPI_File_write(fh, &x[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
  if(mpi->mpiid / mpi->npx == 0)
    MPI_File_write(fh, &xh[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += itot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  if(mpi->mpiid % mpi->npx == 0)
    MPI_File_write(fh, &y[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_Barrier(mpi->commxy);
  fileoff += jtot*sizeof(double);

  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
  if(mpi->mpiid % mpi->npx == 0)
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
  if(mpi->mpiid == 0)
    std::printf("Loading \"%s\"\n", filename);

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

  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
  {
    if(mpi->mpiid == 0)
      std::printf("ERROR \"%s\" cannot be loaded\n", filename);
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

  //if(mpi->mpiid == 0)
  //{
  FILE *pFile;
  pFile = fopen(filename, "rb");
  int n = (2*itot+2*jtot)*sizeof(double);
  fseek(pFile, n, SEEK_SET);
  fread(&z [kstart], sizeof(double), kmax, pFile);
  fread(&zh[kstart], sizeof(double), kmax, pFile);
  fclose(pFile);
  // }

  // calculate the ghost cells
  calculate();

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

  int reqid = 0;

  // first, send the sends
  MPI_Isend(&data[eastout], ncount, eastwestedge, mpi->neast, 1, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;
  MPI_Isend(&data[westout], ncount, eastwestedge, mpi->nwest, 2, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;
  MPI_Isend(&data[northout], ncount, northsouthedge, mpi->nnorth, 1, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;
  MPI_Isend(&data[southout], ncount, northsouthedge, mpi->nsouth, 2, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;

  // second, send the receives, to avoid unnecessary waiting of data
  MPI_Irecv(&data[westin], ncount, eastwestedge, mpi->nwest, 1, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;
  MPI_Irecv(&data[eastin], ncount, eastwestedge, mpi->neast, 2, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;
  MPI_Irecv(&data[southin], ncount, northsouthedge, mpi->nsouth, 1, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;
  MPI_Irecv(&data[northin], ncount, northsouthedge, mpi->nnorth, 2, mpi->commxy, &mpi->reqs[reqid]);
  reqid++;

  MPI_Waitall(reqid, mpi->reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cgrid::transposezx(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;

  int jj = imax;
  int kk = imax*jmax;

  //int kblock = kblock;

  int reqid = 0;

  for(int k=0; k<mpi->npx; k++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + k;

    // determine where to fetch the send data
    int ijks = k*kblock*kk;

    // send the block, tag it with the height (in kblocks) where it should come
    int sendtag = k;
    MPI_Isend(&as[ijks], ncount, transposez, nblock, sendtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  for(int k=0; k<mpi->npx; k++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + k;

    // determine where to store the receive data
    int ijkr = k*jj;

    // and determine what has to be delivered at height k (in kblocks)
    int recvtag = mpi->mpiid % mpi->npx;
    MPI_Irecv(&ar[ijkr], ncount, transposex, nblock, recvtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, mpi->reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cgrid::transposexz(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;

  int jj = imax;
  int kk = imax*jmax;

  // int kblock = kblock;

  int reqid = 0;

  for(int i=0; i<mpi->npx; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + i;

    // determine where to fetch the send data
    int ijks = i*jj;

    // send the block, tag it with the height (in kblocks) where it should come
    int sendtag = mpi->mpiid % mpi->npx;
    MPI_Isend(&as[ijks], ncount, transposex, nblock, sendtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  for(int i=0; i<mpi->npx; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + i;

    // determine where to store the receive data
    int ijkr = i*kblock*kk;

    // and determine what has to be delivered at height i (in kblocks)
    int recvtag = i;
    MPI_Irecv(&ar[ijkr], ncount, transposez, nblock, recvtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }
  MPI_Waitall(reqid, mpi->reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cgrid::transposexy(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;

  int jj = iblock;
  int kk = iblock*jmax;

  int reqid = 0;

  for(int i=0; i<mpi->npy; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid % mpi->npx + i * mpi->npx;

    // determine where to fetch the send data
    int ijks = i*jj;

    // send the block, tag it with the east west location
    int sendtag = i;
    MPI_Isend(&as[ijks], ncount, transposex2, nblock, sendtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  for(int i=0; i<mpi->npy; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid % mpi->npx + i * mpi->npx;

    // determine where to store the receive data
    int ijkr = i*kk;

    // and determine what has to be delivered at depth i (in iblocks)
    int recvtag = mpi->mpiid / mpi->npx;
    MPI_Irecv(&ar[ijkr], ncount, transposey, nblock, recvtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, mpi->reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cgrid::transposeyx(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;

  int jj = iblock;
  int kk = iblock*jmax;

  int reqid = 0;

  for(int i=0; i<mpi->npy; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid % mpi->npx + i * mpi->npx;

    // determine where to fetch the send data
    int ijks = i*kk;

    // send the block, tag it with the east west location
    int sendtag = mpi->mpiid / mpi->npx;
    MPI_Isend(&as[ijks], ncount, transposey, nblock, sendtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  for(int i=0; i<mpi->npy; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid % mpi->npx + i * mpi->npx;

    // determine where to store the receive data
    int ijkr = i*jj;

    // and determine what has to be delivered at depth i (in iblocks)
    int recvtag = i;
    MPI_Irecv(&ar[ijkr], ncount, transposex2, nblock, recvtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, mpi->reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cgrid::transposeyz(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;

  int jj = iblock;
  int kk = iblock*jblock;

  int reqid = 0;

  for(int i=0; i<mpi->npx; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + i;

    // determine where to fetch the send data
    int ijks = i*jblock*jj;

    // send the block, tag it with the height (in kblocks) where it should come
    int sendtag = mpi->mpiid % mpi->npx;
    MPI_Isend(&as[ijks], ncount, transposey2, nblock, sendtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  for(int i=0; i<mpi->npx; i++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + i;

    // determine where to store the receive data
    int ijkr = i*kblock*kk;

    // and determine what has to be delivered at height i (in kblocks)
    int recvtag = i;
    MPI_Irecv(&ar[ijkr], ncount, transposez2, nblock, recvtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }
 
  MPI_Waitall(reqid, mpi->reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cgrid::transposezy(double * restrict ar, double * restrict as)
{
  int nblock;
  int ncount = 1;

  int jj = iblock;
  int kk = iblock*jblock;

  int reqid = 0;

  for(int k=0; k<mpi->npx; k++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + k;

    // determine where to fetch the send data
    int ijks = k*kblock*kk;

    // determine what has to be sent from height i (in kblocks)
    int sendtag = mpi->mpiid % mpi->npx;
    MPI_Isend(&as[ijks], ncount, transposez2, nblock, sendtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  for(int k=0; k<mpi->npx; k++)
  {
    // determine where to send it to
    nblock = mpi->mpiid - mpi->mpiid % mpi->npx + k;

    // determine where to store the receive data
    int ijkr = k*jblock*jj;

    // recv the block, tag is the height (in kblocks) where it should come
    int recvtag = k;
    MPI_Irecv(&ar[ijkr], ncount, transposey2, nblock, recvtag, MPI_COMM_WORLD, &mpi->reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, mpi->reqs, MPI_STATUSES_IGNORE);

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

int cgrid::savefield3d(double * restrict data, char *filename)
{
  MPI_File fh;
  if(MPI_File_open(mpi->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
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

  double *buffer;
  buffer = new double[count];

  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
#pragma ivdep
      for(int i=0; i<imax; i++)
      {
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        ijkb = i + j*jjb + k*kkb;
        buffer[ijkb] = data[ijk];
      }

  fileoff = 0;
  MPI_File_write_at_all(fh, fileoff, buffer, count, MPI_DOUBLE, MPI_STATUS_IGNORE);

  if(MPI_File_close(&fh))
    return 1;

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
  double *buffer;
  buffer = new double[count];

  fileoff = 0;
  MPI_File_read_at_all(fh, fileoff, buffer, count, MPI_DOUBLE, MPI_STATUS_IGNORE);

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
    return 1;

  delete[] buffer;

  return 0;
}
