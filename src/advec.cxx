#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "advec.h"

cadvec::cadvec(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object advec\n");
  grid   = gridin;
  fields = fieldsin;
}

cadvec::~cadvec()
{
  std::printf("Destroying instance of object advec\n");
}

int cadvec::exec()
{
  // advect the flow
  advecu_2nd((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data);
  return 0;
}

// high performance routine
int cadvec::advecu_2nd(double *ut, double *u, double *v, double *w)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dx, dy;
  double *dz;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dx = grid->dx;
  dy = grid->dy;
  dz = grid->dz;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        ut[ijk] += 
              - (  interp2(u[ijk   ] , u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
                 - interp2(u[ijk-ii] , u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) / dx

              - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk-jj])
                 - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) / dy

              - (  interp2(w[ijk+kk-ii], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
                 - interp2(w[ijk   -ii], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) / dz[k];
        //if(i==grid->istart && j == grid->jstart)
        //  std::printf("%3d, %8.5f\n", k, u[k*grid->icells*grid->jcells + j*grid->icells + i]);
      }

  return 0;
}

inline double cadvec::interp2(double a, double b)
{
  double c = 0.5*(a + b);
  return c;
}
