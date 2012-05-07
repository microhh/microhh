#include <cstdio>
#include <mpi.h>
#include "grid.h"
#include "defines.h"
#include "mpiinterface.h"

cmpi::cmpi(cgrid *gridin)
{
  std::printf("Creating instance of object mpi\n");
  grid = gridin;

  initialized = false;
  allocated   = false;
}

cmpi::~cmpi()
{
  if(initialized)
    MPI_Finalize();

  if(allocated)
    delete[] reqs;

  std::printf("Destroying instance of object mpi\n");
}

int cmpi::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&npx, "grid", "npx", 1);
  n += inputin->getItem(&npy, "grid", "npy", 1);

  if(n > 0)
    return 1;
  
  return 0;
}

int cmpi::init()
{
  MPI_Init(NULL, NULL);
  initialized = true;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(nprocs != npx*npy)
  {
    std::printf("ERROR nprocs = %d does not equal npx*npy = %d*%d\n", nprocs, npx, npy);
    return 1;
  }

  int dims    [2] = {npy, npx};
  int periodic[2] = {true, true};

  int n;

  if(MPI_Dims_create(nprocs, 2, dims))
    return 1;

  if(MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, true, &commxy))
    return 1;
  if(MPI_Comm_rank(commxy, &mpiid))
    return 1;

  int mpicoords[2];
  n = MPI_Cart_coords(commxy, mpiid, 2, mpicoords);

  mpicoordx = mpicoords[1];
  mpicoordy = mpicoords[0];

  int dimx[2] = {false, true };
  int dimy[2] = {true , false};

  // MPI_Cart_sub(commxy, dimx, &commx);
  // MPI_Cart_sub(commxy, dimy, &commy);

  if(MPI_Cart_shift(commxy, 1, 1, &nwest , &neast ))
    return 1;
  if(MPI_Cart_shift(commxy, 0, 1, &nsouth, &nnorth))
    return 1;

  // create the MPI types for the cyclic boundary conditions
  int datacount, datablock, datastride;

  // east west
  datacount  = grid->jcells*grid->kcells;
  datablock  = grid->igc;
  datastride = grid->icells;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &eastwestedge);
  MPI_Type_commit(&eastwestedge);

  // north south
  datacount  = grid->kcells;
  datablock  = grid->icells*grid->jgc;
  datastride = grid->icells*grid->jcells;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &northsouthedge);
  MPI_Type_commit(&northsouthedge);

  // transposez
  datacount = grid->imax*grid->jmax*grid->kblock;
  MPI_Type_contiguous(datacount, MPI_DOUBLE, &transposez);
  MPI_Type_commit(&transposez);

  // transposex imax
  datacount  = grid->jmax*grid->kblock;
  datablock  = grid->imax;
  datastride = grid->itot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposex);
  MPI_Type_commit(&transposex);

  // transposex iblock
  datacount  = grid->jmax*grid->kblock;
  datablock  = grid->iblock;
  datastride = grid->itot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposex2);
  MPI_Type_commit(&transposex2);

  // transposey
  datacount  = grid->kblock;
  datablock  = grid->iblock*grid->jmax;
  datastride = grid->iblock*grid->jtot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposey);
  MPI_Type_commit(&transposey);

  // file saving and loading, take C-ordering into account
  int totsize [3] = {grid->kmax, grid->jtot, grid->itot};
  int subsize [3] = {grid->kmax, grid->jmax, grid->imax};
  int substart[3] = {0, mpicoordy*grid->jmax, mpicoordx*grid->imax};
  MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &subarray);
  MPI_Type_commit(&subarray);

  // create the requests arrays for the nonblocking sends
  int npmax;
  npmax = std::max(npx, npy);
  npmax = std::max(npmax, 8);
  reqs = new MPI_Request[npmax*2];

  allocated = true;

  return 0;
} 

int cmpi::boundary_cyclic(double * restrict data)
{
  int ncount = 1;

  // communicate east-west edges
  int eastout = grid->iend-grid->igc;
  int westin  = 0;
  int westout = grid->istart;
  int eastin  = grid->iend;

  int reqid = 0;
  MPI_Isend(&data[eastout], ncount, eastwestedge, neast, 1, commxy, &reqs[reqid]);
  reqid++;
  MPI_Irecv(&data[westin ], ncount, eastwestedge, nwest, 1, commxy, &reqs[reqid]);
  reqid++;
               
  MPI_Isend(&data[westout], ncount, eastwestedge, nwest, 2, commxy, &reqs[reqid]);
  reqid++;
  MPI_Irecv(&data[eastin ], ncount, eastwestedge, neast, 2, commxy, &reqs[reqid]);
  reqid++;

  // communicate north-south edges
  int northout = (grid->jend-grid->jgc)*grid->icells;
  int southin  = 0;
  int southout = grid->jstart*grid->icells;
  int northin  = grid->jend  *grid->icells;

  MPI_Isend(&data[northout], ncount, northsouthedge, nnorth, 1, commxy, &reqs[reqid]);
  reqid++;
  MPI_Irecv(&data[southin ], ncount, northsouthedge, nsouth, 1, commxy, &reqs[reqid]);
  reqid++;

  MPI_Isend(&data[southout], ncount, northsouthedge, nsouth, 2, commxy, &reqs[reqid]);
  reqid++;
  MPI_Irecv(&data[northin ], ncount, northsouthedge, nnorth, 2, commxy, &reqs[reqid]);
  reqid++;

  MPI_Waitall(reqid, reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cmpi::transposezx(double * restrict ar, double * restrict as)
{
  int startk;
  int nblock;
  int ncount = 1;

  int jj = grid->imax;
  int kk = grid->imax*grid->jmax;

  int kblock = grid->kblock;

  int reqid = 0;

  for(int k=0; k<npx; k++)
  {
    // determine where to send it to
    nblock = mpiid - mpiid % npx + k;

    // determine where to fetch the send data
    int ijks = k*kblock*kk;

    // determine where to store the receive data
    int ijkr = k*jj;

    // send the block, tag it with the height (in kblocks) where it should come
    int sendtag = k;
    MPI_Isend(&as[ijks], ncount, transposez, nblock, sendtag, commxy, &reqs[reqid]);
    reqid++;

    // and determine what has to be delivered at height k (in kblocks)
    int recvtag = mpiid % npx;
    MPI_Irecv(&ar[ijkr], ncount, transposex, nblock, recvtag, commxy, &reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cmpi::transposexz(double * restrict ar, double * restrict as)
{
  int startk;
  int nblock;
  int ncount = 1;

  int jj = grid->imax;
  int kk = grid->imax*grid->jmax;

  int kblock = grid->kblock;

  int reqid = 0;

  for(int i=0; i<npx; i++)
  {
    // determine where to send it to
    nblock = mpiid - mpiid % npx + i;

    // determine where to fetch the send data
    int ijks = i*jj;

    // determine where to store the receive data
    int ijkr = i*kblock*kk;

    // send the block, tag it with the height (in kblocks) where it should come
    int sendtag = mpiid % npx;
    MPI_Isend(&as[ijks], ncount, transposex, nblock, sendtag, commxy, &reqs[reqid]);
    reqid++;

    // and determine what has to be delivered at height i (in kblocks)
    int recvtag = i;
    MPI_Irecv(&ar[ijkr], ncount, transposez, nblock, recvtag, commxy, &reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cmpi::transposexy(double * restrict ar, double * restrict as)
{
  int startk;
  int nblock;
  int ncount = 1;

  int jj = grid->iblock;
  int kk = grid->iblock*grid->jmax;

  int reqid = 0;

  for(int i=0; i<npy; i++)
  {
    // determine where to send it to
    nblock = mpiid % npx + i * npx;

    // determine where to fetch the send data
    int ijks = i*jj;

    // determine where to store the receive data
    int ijkr = i*kk;

    // send the block, tag it with the east west location
    int sendtag = i;
    MPI_Isend(&as[ijks], ncount, transposex2, nblock, sendtag, commxy, &reqs[reqid]);
    reqid++;

    // and determine what has to be delivered at depth i (in iblocks)
    int recvtag = mpiid / npx;
    MPI_Irecv(&ar[ijkr], ncount, transposey, nblock, recvtag, commxy, &reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, reqs, MPI_STATUSES_IGNORE);

  return 0;
}

int cmpi::transposeyx(double * restrict ar, double * restrict as)
{
  int startk;
  int nblock;
  int ncount = 1;

  int jj = grid->iblock;
  int kk = grid->iblock*grid->jmax;

  int reqid = 0;

  for(int i=0; i<npy; i++)
  {
    // determine where to send it to
    nblock = mpiid % npx + i * npx;

    // determine where to fetch the send data
    int ijks = i*kk;

    // determine where to store the receive data
    int ijkr = i*jj;

    // send the block, tag it with the east west location
    int sendtag = mpiid / npx;
    MPI_Isend(&as[ijks], ncount, transposey, nblock, sendtag, commxy, &reqs[reqid]);
    reqid++;

    // and determine what has to be delivered at depth i (in iblocks)
    int recvtag = i;
    MPI_Irecv(&ar[ijkr], ncount, transposex2, nblock, recvtag, commxy, &reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, reqs, MPI_STATUSES_IGNORE);
  return 0;
}

int cmpi::transposeyz(double * restrict ar, double * restrict as)
{
  int startk;
  int nblock;
  int ncount = 1;

  int jj = grid->iblock;
  int kk = grid->iblock*grid->jmax;

  int reqid = 0;

  for(int i=0; i<npy; i++)
  {
    // determine where to send it to
    nblock = mpiid % npx + i * npx;

    // determine where to fetch the send data
    int ijks = i*kk;

    // determine where to store the receive data
    int ijkr = i*jj;

    // send the block, tag it with the east west location
    int sendtag = mpiid / npx;
    MPI_Isend(&as[ijks], ncount, transposey, nblock, sendtag, commxy, &reqs[reqid]);
    reqid++;

    // and determine what has to be delivered at depth i (in iblocks)
    int recvtag = i;
    MPI_Irecv(&ar[ijkr], ncount, transposez, nblock, recvtag, commxy, &reqs[reqid]);
    reqid++;
  }

  MPI_Waitall(reqid, reqs, MPI_STATUSES_IGNORE);
  return 0;
}

int cmpi::getmax(double *var)
{
  double varl = *var;
  MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_MAX, commxy);

  return 0;
}

int cmpi::getsum(double *var)
{
  double varl = *var;
  MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_SUM, commxy);

  return 0;
}

int cmpi::writefield3d(double *data, char *filename)
{
  int n = 0;
  MPI_File fh;
  if(MPI_File_open(commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, NULL, MPI_INFO_NULL);

  // extract the data from the 3d field without the ghost cells
  int ijk, jj, kk;
  jj  = grid->icells;
  kk  = grid->icells*grid->jcells;
  ijk = grid->istart + grid->jstart*jj + grid->kstart*kk;

  int count = grid->imax * grid->jmax * grid->kmax;
  for(int n=0; n<count; n++)
    data[n] = mpiid;

  fileoff = mpicoordx*grid->imax + mpicoordy*grid->itot*grid->jmax;
  MPI_File_write_at_all(fh, fileoff, data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);

  if(MPI_File_close(&fh))
    return 1;

  return 0;
}

int cmpi::readfield3d(double *data, char *filename)
{
  return 0;
}
