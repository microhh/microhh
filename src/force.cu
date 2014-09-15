#include "force.h"
#include "grid.h"
#include "fields.h"

__global__ void force_flux_step1(double * const __restrict__ usum, double * const __restrict__ utsum,
                                 const double * const __restrict__ u, const double * const __restrict__ ut,
                                 const double * const __restrict__ dz,
                                 const int jj, const int kk, 
                                 const int istart, const int jstart, const int kstart,
                                 const int iend,   const int jend,   const int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    usum [ijk] = u [ijk]*dz[k];
    utsum[ijk] = ut[ijk]*dz[k];
  }
}

__global__ void force_flux_step2(double * const __restrict__ ut,
                                 const double fbody,
                                 const int jj, const int kk, 
                                 const int istart, const int jstart, const int kstart,
                                 const int iend,   const int jend,   const int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] += fbody;
  }
}

#ifdef USECUDA
int cforce::exec(double dt)
{
  if(swlspres == "uflux")
  {
    const int blocki = 128;
    const int blockj = 2;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    force_flux_step1<<<gridGPU, blockGPU>>>(&fields->a["tmp1"]->data_g[offs], &fields->a["tmp2"]->data_g[offs],
                                            &fields->u->data_g[offs], &fields->ut->data_g[offs],
                                            grid->dz_g,
                                            grid->icellsp, grid->ijcellsp,
                                            grid->istart,  grid->jstart, grid->kstart,
                                            grid->iend,    grid->jend,   grid->kend);

    double uavg  = grid->getsum_g(&fields->a["tmp1"]->data_g[offs], fields->a["tmp3"]->data_g); 
    double utavg = grid->getsum_g(&fields->a["tmp2"]->data_g[offs], fields->a["tmp3"]->data_g); 

    uavg  = uavg  / (grid->itot*grid->jtot*grid->zsize);
    utavg = utavg / (grid->itot*grid->jtot*grid->zsize);

    double fbody = (uflux - uavg - grid->utrans) / dt - utavg;

    force_flux_step2<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs],
                                            fbody,
                                            grid->icellsp, grid->ijcellsp,
                                            grid->istart,  grid->jstart, grid->kstart,
                                            grid->iend,    grid->jend,   grid->kend);
  }

  /*
  else if(swlspres == "geo")
  {
    if(grid->swspatialorder == "2")
      coriolis_2nd(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
    else if(grid->swspatialorder == "4")
      coriolis_4th(fields->ut->data, fields->vt->data, fields->u->data, fields->v->data, ug, vg);
  }

  if(swls == "1")
  {
    for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
      lssource(fields->st[*it]->data, lsprofs[*it]);
  }

  if(swwls == "1")
  {
    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      advecwls_2nd(it->second->data, fields->s[it->first]->datamean, wls, grid->dzhi);
  }
  */

  return 0;
}
#endif
