#include "advec_2.h"
#include "grid.h"
#include "fields.h"
#include "tools.h"

__device__ double interp2(double a, double b)
{
  return 0.5*(a + b);
}

__global__ void advec_2_advecu(double * __restrict__ ut, double * __restrict__ u, 
                               double * __restrict__ v, double * __restrict__ w,
                               double * __restrict__ dzi, double dxi, double dyi, 
                               int jj, int kk, int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] += 
          - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
             - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

          - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
             - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi 

          - (  interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
             - interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) * dzi[k];
  }
}

__global__ void advec_2_advecv(double * __restrict__ vt, double * __restrict__ u, 
                               double * __restrict__ v, double * __restrict__ w,
                               double * __restrict__ dzi, double dxi, double dyi, 
                               int jj, int kk, int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    vt[ijk] += 
          - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
             - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

          - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
             - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

          - (  interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
             - interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) * dzi[k];
  }
}

__global__ void advec_2_advecw(double * __restrict__ wt, double * __restrict__ u, 
                               double * __restrict__ v, double * __restrict__ w,
                               double * __restrict__ dzhi, double dxi, double dyi, 
                               int jj, int kk, int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart + 1;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    wt[ijk] += 
          - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
             - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

          - (  interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
             - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

          - (  interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
             - interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) * dzhi[k];
  }
}

__global__ void advec_2_advecuvw(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                                 double * __restrict__ u,  double * __restrict__ v,  double * __restrict__ w,
                                 double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi, 
                                 int jj, int kk, int istart, int jstart, int kstart,
                                 int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] += 
          - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
             - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

          - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
             - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi 

          - (  interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
             - interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) * dzi[k];

    vt[ijk] += 
          - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
             - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

          - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
             - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

          - (  interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
             - interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) * dzi[k];

    if(k>kstart)
    {
      wt[ijk] += 
            - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
               - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

            - (  interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
               - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

            - (  interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
               - interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) * dzhi[k];
    }
  }
}

// !!!!!!!!!! at the moment (for testing) this only works if itot & jtot are multiple of blockdims !!!!!!!! 
__global__ void advec_2_advecuvw_smem(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                                      double * __restrict__ u,  double * __restrict__ v,  double * __restrict__ w,
                                      double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi, 
                                      int jj, int kk, int istart, int jstart, int kstart,
                                      int iend, int jend, int kend, int igc, int jgc)
{
  extern __shared__ double shared[];

  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    const int of = (blockDim.x + 2 * igc) * (blockDim.y + 2 * jgc); // Stride per variable in __shared__
    const int ii = 1; // Stride distance x in __global__ and __shared__
    const int jjs = blockDim.x + 2 * igc; // Stride distance y in __shared__
    const int ijk = i + j*jj + k*kk; // Index in __global__ array
    const int ijs = threadIdx.x + igc + (blockDim.x + 2 * igc) * (threadIdx.y + jgc); // Index in __shared__ array

    double* us = (double*)shared; 
    double* vs = (double*)&shared[1*of];
    double* ws = (double*)&shared[2*of];

    // Each thread read variable into shared array
    us[ijs] = u[ijk];
    vs[ijs] = v[ijk];
    ws[ijs] = w[ijk];

    // Some threads read ghost cells (how efficent is this?)
    if(threadIdx.x == 0 && threadIdx.y == 0) {
      us[ijs-ii-jjs] = u[ijk-ii-jj]; 
      vs[ijs-ii-jjs] = v[ijk-ii-jj]; 
      ws[ijs-ii-jjs] = w[ijk-ii-jj]; 
    }
    if(threadIdx.x == 0 && threadIdx.y == blockDim.y-1) {
      us[ijs-ii+jjs] = u[ijk-ii+jj]; 
      vs[ijs-ii+jjs] = v[ijk-ii+jj]; 
      ws[ijs-ii+jjs] = w[ijk-ii+jj]; 
    }
    if(threadIdx.x == blockDim.x-1 && threadIdx.y == 0) {
      us[ijs+ii-jjs] = u[ijk+ii-jj]; 
      vs[ijs+ii-jjs] = v[ijk+ii-jj]; 
      ws[ijs+ii-jjs] = w[ijk+ii-jj]; 
    }
    if(threadIdx.x == blockDim.x-1 && threadIdx.y == blockDim.y-1) {
      us[ijs+ii+jjs] = u[ijk+ii+jj]; 
      vs[ijs+ii+jjs] = v[ijk+ii+jj]; 
      ws[ijs+ii+jjs] = w[ijk+ii+jj]; 
    }
    if(threadIdx.x == 0) {
      us[ijs-ii]  = u[ijk-ii];
      vs[ijs-ii]  = v[ijk-ii];
      ws[ijs-ii]  = w[ijk-ii];
    } 
    if(threadIdx.x == blockDim.x-1) {
      us[ijs+ii]  = u[ijk+ii]; 
      vs[ijs+ii]  = v[ijk+ii]; 
      ws[ijs+ii]  = w[ijk+ii]; 
    }
    if(threadIdx.y == 0) {
      us[ijs-jjs] = u[ijk-jj]; 
      vs[ijs-jjs] = v[ijk-jj]; 
      ws[ijs-jjs] = w[ijk-jj]; 
    }
    if(threadIdx.y == blockDim.y-1) {
      us[ijs+jjs] = u[ijk+jj]; 
      vs[ijs+jjs] = v[ijk+jj]; 
      ws[ijs+jjs] = w[ijk+jj]; 
    }

    __syncthreads(); 

    ut[ijk] += 
          - (  interp2(us[ijs       ], us[ijs+ii ]) * interp2(us[ijs    ], us[ijs+ii ])
             - interp2(us[ijs-ii    ], us[ijs    ]) * interp2(us[ijs-ii ], us[ijs    ]) ) * dxi

          - (  interp2(vs[ijs-ii+jjs], vs[ijs+jjs]) * interp2(us[ijs    ], us[ijs+jjs])
             - interp2(vs[ijs-ii    ], vs[ijs    ]) * interp2(us[ijs-jjs], us[ijs    ]) ) * dyi 

          - (  interp2(w[ijk-ii+kk  ], w[ijk+kk  ]) * interp2(us[ijs    ], u[ijk+kk  ])
             - interp2(ws[ijs-ii    ], ws[ijs    ]) * interp2(u[ijk-kk  ], us[ijs    ]) ) * dzi[k];

    vt[ijk] += 
          - (  interp2(us[ijs+ii-jjs], us[ijs+ii ]) * interp2(vs[ijs    ], vs[ijs+ii ])
             - interp2(us[ijs-jjs   ], us[ijs    ]) * interp2(vs[ijs-ii ], vs[ijs    ]) ) * dxi

          - (  interp2(vs[ijs       ], vs[ijs+jjs]) * interp2(vs[ijs    ], vs[ijs+jjs])
             - interp2(vs[ijs-jjs   ], vs[ijs    ]) * interp2(vs[ijs-jjs], vs[ijs    ]) ) * dyi

          - (  interp2(w[ijk-jj+kk  ], w[ijk+kk  ]) * interp2(vs[ijs    ], v[ijk+kk  ])
             - interp2(ws[ijs-jjs   ], ws[ijs    ]) * interp2(v[ijk-kk  ], vs[ijs    ]) ) * dzi[k];

    if(k>kstart)
    {
      wt[ijk] += 
            - (  interp2(u[ijk+ii-kk], us[ijs+ii ]) * interp2(ws[ijs    ], ws[ijs+ii ])
               - interp2(u[ijk   -kk], us[ijs    ]) * interp2(ws[ijs-ii ], ws[ijs    ]) ) * dxi

            - (  interp2(v[ijk+jj-kk], vs[ijs+jjs]) * interp2(ws[ijs    ], ws[ijs+jjs])
               - interp2(v[ijk   -kk], vs[ijs    ]) * interp2(ws[ijs-jjs], ws[ijs    ]) ) * dyi

            - (  interp2(ws[ijs     ], w[ijk+kk  ]) * interp2(ws[ijs    ], w[ijk+kk  ])
               - interp2(w[ijk-kk   ], ws[ijs    ]) * interp2(w[ijk-kk  ], ws[ijs    ]) ) * dzhi[k];
    }
  }
}

__global__ void advec_2_advecs(double * __restrict__ st, double * __restrict__ s, 
                               double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                               double * __restrict__ dzi, double dxi, double dyi, 
                               int jj, int kk, int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    st[ijk] += 
          - (  u[ijk+ii] * interp2(s[ijk   ], s[ijk+ii])
             - u[ijk   ] * interp2(s[ijk-ii], s[ijk   ]) ) * dxi

          - (  v[ijk+jj] * interp2(s[ijk   ], s[ijk+jj])
             - v[ijk   ] * interp2(s[ijk-jj], s[ijk   ]) ) * dyi 

          - (  w[ijk+kk] * interp2(s[ijk   ], s[ijk+kk])
             - w[ijk   ] * interp2(s[ijk-kk], s[ijk   ]) ) * dzi[k];
  }
}

__global__ void advec_2_calccfl(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
                                double * __restrict__ tmp1, double * __restrict__ dzi, double dxi, double dyi, 
                                int jj, int kk, int istart, int jstart, int kstart,
                                int iend, int jend, int kend,
                                int icells, int jcells, int kcells)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 
  int k = blockIdx.z; 
  int ii = 1;
  int ijk = i + j*jj + k*kk;

  if(i >= istart && i < iend && j >= jstart && j < jend && k >= kstart && k < kend)
    tmp1[ijk] = std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + 
                std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + 
                std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k];
  else if(i < icells && j < jcells && k < kcells) 
    tmp1[ijk] = 0.;
}

#ifdef USECUDA
int cadvec_2::exec()
{
  //fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;
  
  const int nsmem = (blocki + 2 * grid->igc) * (blockj + 2 * grid->jgc);  

  //advec_2_advecu<<<gridGPU, blockGPU>>>(fields->ut->data_g, fields->u->data_g, fields->v->data_g, 
  //                                      fields->w->data_g, grid->dzi_g, dxi, dyi,
  //                                      grid->icells, grid->ijcells,
  //                                      grid->istart, grid->jstart, grid->kstart,
  //                                      grid->iend,   grid->jend, grid->kend);

  //advec_2_advecv<<<gridGPU, blockGPU>>>(fields->vt->data_g, fields->u->data_g, fields->v->data_g, 
  //                                      fields->w->data_g, grid->dzi_g, dxi, dyi,
  //                                      grid->icells, grid->ijcells,
  //                                      grid->istart, grid->jstart, grid->kstart,
  //                                      grid->iend,   grid->jend, grid->kend);

  //advec_2_advecw<<<gridGPU, blockGPU>>>(fields->wt->data_g, fields->u->data_g, fields->v->data_g, 
  //                                      fields->w->data_g, grid->dzhi_g, dxi, dyi,
  //                                      grid->icells, grid->ijcells,
  //                                      grid->istart, grid->jstart, grid->kstart,
  //                                      grid->iend,   grid->jend, grid->kend);

  advec_2_advecuvw<<<gridGPU, blockGPU>>>(fields->ut->data_g, fields->vt->data_g, fields->wt->data_g, 
                                          fields->u->data_g,  fields->v->data_g,  fields->w->data_g, 
                                          grid->dzi_g, grid->dzhi_g, dxi, dyi,
                                          grid->icells, grid->ijcells,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend, grid->kend);

  //advec_2_advecuvw_smem<<<gridGPU, blockGPU, 3*nsmem*sizeof(double)>>>
  //                                        (fields->ut->data_g, fields->vt->data_g, fields->wt->data_g, 
  //                                         fields->u->data_g,  fields->v->data_g,  fields->w->data_g, 
  //                                         grid->dzi_g, grid->dzhi_g, dxi, dyi,
  //                                         grid->icells, grid->ijcells,
  //                                         grid->istart, grid->jstart, grid->kstart,
  //                                         grid->iend,   grid->jend, grid->kend,
  //                                         grid->igc,    grid->jgc);

  for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advec_2_advecs<<<gridGPU, blockGPU>>>((*it->second).data_g, (*fields->s[it->first]).data_g, 
                                          fields->u->data_g, fields->v->data_g, fields->w->data_g, 
                                          grid->dzi_g, dxi, dyi,
                                          grid->icells, grid->ijcells,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend, grid->kend);

  //fields->backwardGPU();
  return 0;
}
#endif

#ifdef USECUDA
double cadvec_2::calccfl(double * u, double * v, double * w, double * dzi, double dt)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);
  double cfl = 0;

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  //fields->forwardGPU();

  advec_2_calccfl<<<gridGPU, blockGPU>>>(fields->u->data_g, fields->v->data_g, fields->w->data_g, 
                                         fields->a["tmp1"]->data_g, grid->dzi_g, dxi, dyi,
                                         grid->icells, grid->ijcells,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->iend,   grid->jend, grid->kend,
                                         grid->icells, grid->jcells, grid->kcells);

  cfl = maximum_gpu(fields->a["tmp1"]->data_g, grid->ncells);
  grid->getmax(&cfl);
  cfl = cfl*dt;

  return cfl;
}
#endif


