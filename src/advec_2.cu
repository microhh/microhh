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
  
  const int offs = grid->memoffset;

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

  advec_2_advecuvw<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs], 
                                          &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs], 
                                          grid->dzi_g, grid->dzhi_g, dxi, dyi,
                                          grid->icellsp, grid->ijcellsp,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend, grid->kend);
  
  for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advec_2_advecs<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->s[it->first]->data_g[offs], 
                                          &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs], 
                                          grid->dzi_g, dxi, dyi,
                                          grid->icellsp, grid->ijcellsp,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend, grid->kend);

  //fields->backwardGPU();

  return 0;
}
#endif


#ifdef USECUDA
double cadvec_2::calccfl(double * u, double * v, double * w, double * dzi, double dt)
{
  //fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icellsp/blocki + (grid->icellsp%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);
  double cfl = 0;

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  const int offs = grid->memoffset;

  /* TODO: write algoritm that reduces only the domain itself, excluding ghost and padding cells
     With padding, the overhead of the ghost/padding cells can get quite large */
  advec_2_calccfl<<<gridGPU, blockGPU>>>(&fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs], 
                                         &fields->a["tmp1"]->data_g[offs], grid->dzi_g, dxi, dyi,
                                         grid->icellsp, grid->ijcellsp,
                                         grid->istart,  grid->jstart, grid->kstart,
                                         grid->iend,    grid->jend,   grid->kend,
                                         grid->icellsp, grid->jcells, grid->kcells);

  cfl = maximum_gpu(&fields->a["tmp1"]->data_g[offs], grid->ncellsp-offs);
  grid->getmax(&cfl);
  cfl = cfl*dt;

  return cfl;
}
#endif


