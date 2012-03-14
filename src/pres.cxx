#include <cstdio>
#include <cmath>
#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"

cpres::cpres(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object pres\n");
  grid   = gridin;
  fields = fieldsin;
}

cpres::~cpres()
{
  std::printf("Destroying instance of object pres\n");
}

int cpres::exec(double dt)
{
  // cyclic boundaries for tendencies 
  (*fields->ut).boundary_cyclic();
  (*fields->vt).boundary_cyclic();
  (*fields->wt).boundary_cyclic();

  // create the input for the pressure solver
  pres_2nd_in((*fields->p ).data,
              (*fields->u ).data, (*fields->v ).data, (*fields->w ).data,
              (*fields->ut).data, (*fields->vt).data, (*fields->wt).data, 
              grid->dzi, dt);

  // solve the system
  pres_2nd_solve((*fields->p).data);
  
  // get the pressure tendencies from the pressure field
  pres_2nd_out((*fields->ut).data, (*fields->vt).data, (*fields->wt).data, 
               (*fields->p ).data, grid->dzhi);

  return 0;
}

int cpres::init()
{
  pres_2nd_init();

  return 0;
}

int cpres::pres_2nd_init()
{
  int itot, jtot;
  itot = grid->itot;
  jtot = grid->jtot;

  fftini  = new double[itot];
  fftouti = new double[itot];
  fftinj  = new double[jtot];
  fftoutj = new double[jtot];

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_PATIENT);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_PATIENT);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_PATIENT);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_PATIENT);

  bmati = new double[itot];
  bmatj = new double[itot];
  
  // compute the modified wave numbers of the 2nd order scheme
  double dxidxi = 1./(grid->dx*grid->dx);
  double dyidyi = 1./(grid->dy*grid->dy);

  const double pi = std::acos(-1.);

  for(int j=0; j<jtot/2+1; j++)
    bmatj[j] = 2. * (std::cos(2.*pi*j/jtot)-1.) * dyidyi;

  for(int j=jtot/2+1; j<jtot; j++)
    bmatj[j] = bmatj[jtot-j];

  for(int i=0; i<itot/2+1; i++)
    bmati[i] = 2. * (std::cos(2.*pi*i/itot)-1.) * dxidxi;

  for(int i=itot/2+1; i<itot; i++)
    bmati[i] = bmati[itot-i];

  // for(int i=0; i<itot; i++)
  //  std::printf("%d, %f\n", i, bmatj[i]);

  return 0;
}

/*  subroutine initpres_2nd
    use modgriddata,  only : kmax, itot, jtot, dxi, dyi, dz, dzhi
    use modconstants, only : pi

    integer :: i, j, k

    ! prepare FFTs
    allocate(fftini(itot))
    allocate(fftouti(itot))
    allocate(fftinj(jtot))
    allocate(fftoutj(jtot))

    iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, fftw_r2hc, fftw_patient)
    iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, fftw_hc2r, fftw_patient)
    jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, fftw_r2hc, fftw_patient)
    jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, fftw_hc2r, fftw_patient)

    ! allocate matrix coefficients
    allocate(bmati(itot))
    allocate(bmatj(jtot))

    allocate(a  (1:kmax+1))
    allocate(b  (0:kmax+1))
    allocate(c  (0:kmax))
    allocate(blp(1:kmax))

    allocate(x    (1:kmax  ))
    allocate(alp  (2:kmax  ))
    allocate(clp  (1:kmax-1))

    a(:) = 0.
    b(:) = 0.
    c(:) = 0.
    x(:) = 0.

    ! Create coefficients for b-vector in tridiagonal solver
    do j = 1, jtot/2+1
      bmatj(j) = 2. * (cos(2. * pi * (j-1) / jtot) - 1.) * dyi * dyi
    end do
    do j = jtot/2+2, jtot
      bmatj(j) = bmatj(jtot - j + 2)
    end do

    do i = 1, itot/2+1
      bmati(i) = 2. * (cos(2. * pi * (i-1) / itot) - 1.) * dxi * dxi
    end do

    do i = itot/2+2, itot
      bmati(i) = bmati(itot - i + 2)
    end do

    ! create vectors that go into the tridiagonal matrix solver
    do k = 1, kmax
      a(k)   = dz(k) * dzhi(k)
      c(k)   = dz(k) * dzhi(k+1)
    end do

    ! fill in boundary conditions in matrix, no pressure gradient on top and
    ! bottom
    c(0)      = dz(1) * dzhi(1)
    a(kmax+1) = -dz(kmax) * dzhi(kmax+1)
  end subroutine initpres_2nd
*/

int cpres::pres_2nd_in(double * __restrict__ p, 
                       double * __restrict__ u , double * __restrict__ v , double * __restrict__ w , 
                       double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                       double * __restrict__ dzi,
                       double dt)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dxi, dyi;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        p[ijk] = ( (ut[ijk+ii] + u[ijk+ii] / dt) - (ut[ijk] + u[ijk] / dt) ) * dxi
               + ( (vt[ijk+jj] + v[ijk+jj] / dt) - (vt[ijk] + v[ijk] / dt) ) * dyi
               + ( (wt[ijk+kk] + w[ijk+kk] / dt) - (wt[ijk] + w[ijk] / dt) ) * dzi[k];
           
      }

  return 0;
}

int cpres::pres_2nd_solve(double * __restrict__ p)
{
  int i,j,k,ii,jj,kk,ijk;
  int imax, jmax, kmax, itot, jtot, ktot;

  imax = grid->imax;
  jmax = grid->jmax;
  kmax = grid->kmax;
  itot = grid->itot;
  jtot = grid->jtot;
  ktot = grid->ktot;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // do the first fourier transform
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
    {
      for(int i=0;i<itot;i++)
      { 
        ijk = i + j*jj + k*kk;
        fftini[i] = p[ijk];
      }

      fftw_execute_r2r(iplanf, fftini, fftouti);

      for(int i=0;i<itot;i++)
      {
        ijk = i + j*jj+ k*kk;
        p[ijk] = fftouti[i];
      }
    }

  // TRANSPOSE

  // do the second fourier transform
  for(int k=0; k<kmax; k++)
    for(int i=0; i<imax; i++)
    {
      for(int j=0;j<jtot;j++)
      { 
        ijk = i + j*jj + k*kk;
        fftinj[j] = p[ijk];
      }

      fftw_execute_r2r(jplanf, fftinj, fftoutj);

      for(int j=0;j<jtot;j++)
      {
        ijk = i + j*jj+ k*kk;
        p[ijk] = fftoutj[j];
      }
    }

  // solve the tridiagonal system

  // TRANSPOSE
  
  // transform the second transform back
  for(int k=0; k<kmax; k++)
    for(int i=0; i<imax; i++)
    {
      for(int j=0;j<jtot;j++)
      { 
        ijk = i + j*jj + k*kk;
        fftinj[j] = p[ijk];
      }

      fftw_execute_r2r(jplanb, fftinj, fftoutj);

      for(int j=0;j<jtot;j++)
      {
        ijk = i + j*jj+ k*kk;
        p[ijk] = fftoutj[j];
      }
    }

    // TRANSPOSE
    
    // transform the first transform back
    for(int k=0; k<kmax; k++)
      for(int j=0; j<jmax; j++)
      {
        for(int i=0;i<itot;i++)
        { 
          ijk = i + j*jj + k*kk;
          fftini[i] = p[ijk];
        }

        fftw_execute_r2r(iplanb, fftini, fftouti);

        for(int i=0;i<itot;i++)
        {
          ijk = i + j*jj+ k*kk;
          p[ijk] = fftouti[i];
        }
      }


  /*!! Solve the poisson equation per wave number
    do j = 1, jmax
      do i = 1, imax ! rest of array is filled with zeros
        iindex = mpicoordx * imax + i
        jindex = mpicoordy * jmax + j

        ! create vectors that go into the tridiagonal matrix solver
        do k = 1, kmax
          b(k) = dz(k) * dz(k) * (bmati(iindex) + bmatj(jindex)) - (a(k) + c(k))
          x(k) = dz(k) * dz(k) * p(k,i,j)
        end do

        ! and now.... solve the matrix!
        blp(:) = b(1:kmax)

        ! substitute BC's
        blp(1)    = blp(1) + a(1)

        if(iindex == 1 .and. jindex == 1) then
          ! for wave number 0, which contains average, set pressure at top to zero
          blp(kmax) = blp(kmax) - c(kmax)
        else
          ! set dp/dz = 0
          blp(kmax) = blp(kmax) + c(kmax)
        end if

        ! LAPACK band storage
        alp(:) = a(2:kmax)
        clp(:) = c(1:kmax-1)

        ! call LAPACK tridiagonal matrix solver
        call dgtsv(kmax, 1, alp, blp, clp, x, kmax, info)
        
        ! update the pressure (in fourier space, still)
        p(1:kmax,i,j) = x(1:kmax)
      end do
    end do

    ! recompute the tendencies
    call boundary_cyclic(p)

    ! no pressure gradient at the surface
    p(0,:,:) = p(1,:,:) */

  return 0;
}

int cpres::pres_2nd_out(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                        double * __restrict__ p , double * __restrict__ dzhi)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dxi, dyi;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        ut[ijk] = ut[ijk] - (p[ijk] - p[ijk-ii]) * dxi;
        vt[ijk] = vt[ijk] - (p[ijk] - p[ijk-jj]) * dyi;
        wt[ijk] = wt[ijk] - (p[ijk] - p[ijk-kk]) * dzhi[k];
      }

  return 0;
}

