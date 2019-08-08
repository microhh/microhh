! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Kernels for computing broadband fluxes by summing over all elements in the spectral dimension
!
! -------------------------------------------------------------------------------------------------
module mo_fluxes_broadband_kernels
  use, intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp
  implicit none
  private
  public :: sum_broadband, net_broadband

  interface net_broadband
    module procedure net_broadband_full, net_broadband_precalc
  end interface net_broadband
contains
  ! ----------------------------------------------------------------------------
    !
    ! Spectral reduction over all points
    !
  pure subroutine sum_broadband(ncol, nlev, ngpt, spectral_flux, broadband_flux) bind(C, name="sum_broadband")
    integer,                               intent(in ) :: ncol, nlev, ngpt
    real(wp), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
    real(wp), dimension(ncol, nlev),       intent(out) :: broadband_flux

    integer :: icol, ilev, igpt

    !$acc enter data copyin(spectral_flux) create(broadband_flux)
    !$acc parallel loop collapse(2)
    do ilev = 1, nlev
      do icol = 1, ncol
        broadband_flux(icol, ilev) =  spectral_flux(icol, ilev, 1)
      end do
    end do
    !$acc parallel loop collapse(3)
    do igpt = 2, ngpt
      do ilev = 1, nlev
        do icol = 1, ncol
          !$acc atomic update
          broadband_flux(icol, ilev) = broadband_flux(icol, ilev) + spectral_flux(icol, ilev, igpt)
        end do
      end do
    end do
    !$acc exit data delete(spectral_flux) copyout(broadband_flux)
  end subroutine sum_broadband
  ! ----------------------------------------------------------------------------
  !
  ! Net flux: Spectral reduction over all points
  !
  pure subroutine net_broadband_full(ncol, nlev, ngpt, spectral_flux_dn, spectral_flux_up, broadband_flux_net) &
    bind(C, name="net_broadband_full")
    integer,                               intent(in ) :: ncol, nlev, ngpt
    real(wp), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux_dn, spectral_flux_up
    real(wp), dimension(ncol, nlev),       intent(out) :: broadband_flux_net

    integer  :: icol, ilev, igpt
    real(wp) :: diff

    !$acc enter data copyin(spectral_flux_dn, spectral_flux_up) create(broadband_flux_net)
    !$acc parallel loop collapse(2)
    do ilev = 1, nlev
      do icol = 1, ncol
        diff = spectral_flux_dn(icol, ilev, 1   ) - spectral_flux_up(icol, ilev,     1)
        broadband_flux_net(icol, ilev) = diff
      end do
    end do
    !$acc parallel loop collapse(3)
    do igpt = 2, ngpt
      do ilev = 1, nlev
        do icol = 1, ncol
          diff = spectral_flux_dn(icol, ilev, igpt) - spectral_flux_up(icol, ilev, igpt)
          !$acc atomic update
          broadband_flux_net(icol, ilev) = broadband_flux_net(icol, ilev) + diff
        end do
      end do
    end do
    !$acc exit data delete(spectral_flux_dn, spectral_flux_up) copyout(broadband_flux_net)
  end subroutine net_broadband_full
  ! ----------------------------------------------------------------------------
  !
  ! Net flux when bradband flux up and down are already available
  !
  pure subroutine net_broadband_precalc(ncol, nlev, flux_dn, flux_up, broadband_flux_net) &
    bind(C, name="net_broadband_precalc")
    integer,                         intent(in ) :: ncol, nlev
    real(wp), dimension(ncol, nlev), intent(in ) :: flux_dn, flux_up
    real(wp), dimension(ncol, nlev), intent(out) :: broadband_flux_net

    integer  :: icol, ilev
    !$acc enter data copyin(flux_dn, flux_up) create(broadband_flux_net)
    !$acc parallel loop collapse(2)
    do ilev = 1, nlev
      do icol = 1, ncol
         broadband_flux_net(icol,ilev) = flux_dn(icol,ilev) - flux_up(icol,ilev)
       end do
    end do
    !$acc exit data delete(flux_dn, flux_up) copyout(broadband_flux_net)
  end subroutine net_broadband_precalc
  ! ----------------------------------------------------------------------------
end module mo_fluxes_broadband_kernels
