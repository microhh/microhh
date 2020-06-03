! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2019,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
module mo_rte_util_array
!
! This module provide utilites for sanitizing input arrays:
!    checking values and sizes
! These are in a module so code can be written for both CPUs and GPUs
! Used only by Fortran classes so routines don't need C bindings and can use assumed-shape
!
  use mo_rte_kind,      only: wp, wl
  implicit none
  interface any_vals_less_than
    module procedure any_vals_less_than_1D,        any_vals_less_than_2D,        any_vals_less_than_3D
    module procedure any_vals_less_than_1D_masked, any_vals_less_than_2D_masked, any_vals_less_than_3D_masked
  end interface
  interface any_vals_outside
    module procedure any_vals_outside_1D,        any_vals_outside_2D,        any_vals_outside_3D
    module procedure any_vals_outside_1D_masked, any_vals_outside_2D_masked, any_vals_outside_3D_masked
  end interface
  interface zero_array
    module procedure zero_array_1D, zero_array_3D, zero_array_4D
  end interface
  interface extents_are
    module procedure extents_are_1D, extents_are_2D, extents_are_3D
    module procedure extents_are_4D, extents_are_5D, extents_are_6D
    module procedure extents_are_2d_int
  end interface extents_are

  private
  public :: any_vals_less_than, any_vals_outside, extents_are, zero_array
contains
  !-------------------------------------------------------------------------------------------------
  ! Values less than a floor
  !-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_1D(array, check_value)
    real(wp), dimension(:), intent(in) :: array
    real(wp),               intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    !$acc end kernels

    any_vals_less_than_1D = (minValue < check_value)

  end function any_vals_less_than_1D
!-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_2D(array, check_value)
    real(wp), dimension(:,:), intent(in) :: array
    real(wp),                 intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    !$acc end kernels

    any_vals_less_than_2D = (minValue < check_value)

  end function any_vals_less_than_2D
!-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_3D(array, check_value)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    !$acc end kernels

    any_vals_less_than_3D = (minValue < check_value)

  end function any_vals_less_than_3D
  !-------------------------------------------------------------------------------------------------
  ! Masked versions
  !-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_1D_masked(array, mask, check_value)
    real   (wp), dimension(:), intent(in) :: array
    logical(wl), dimension(:), intent(in) :: mask
    real   (wp),               intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array, mask=mask)
    !$acc end kernels

    any_vals_less_than_1D_masked = (minValue < check_value)

  end function any_vals_less_than_1D_masked
  !-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_2D_masked(array, mask, check_value)
    real   (wp), dimension(:,:), intent(in) :: array
    logical(wl), dimension(:,:), intent(in) :: mask
    real   (wp),                 intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array, mask=mask)
    !$acc end kernels

    any_vals_less_than_2D_masked = (minValue < check_value)

  end function any_vals_less_than_2D_masked
  !-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_3D_masked(array, mask, check_value)
    real   (wp), dimension(:,:,:), intent(in) :: array
    logical(wl), dimension(:,:,:), intent(in) :: mask
    real   (wp),                   intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array, mask=mask)
    !$acc end kernels

    any_vals_less_than_3D_masked = (minValue < check_value)

  end function any_vals_less_than_3D_masked
  !-------------------------------------------------------------------------------------------------
  ! Values outside a range
  !-------------------------------------------------------------------------------------------------
  logical function any_vals_outside_1D(array, checkMin, checkMax)
    real(wp), dimension(:), intent(in) :: array
    real(wp),               intent(in) :: checkMin, checkMax

    real(wp) :: minValue, maxValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    maxValue = maxval(array)
    !$acc end kernels
    any_vals_outside_1D = minValue < checkMin .or. maxValue > checkMax

  end function any_vals_outside_1D
! ----------------------------------------------------------
  logical function any_vals_outside_2D(array, checkMin, checkMax)
    real(wp), dimension(:,:), intent(in) :: array
    real(wp),                 intent(in) :: checkMin, checkMax

    real(wp) :: minValue, maxValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    maxValue = maxval(array)
    !$acc end kernels
    any_vals_outside_2D = minValue < checkMin .or. maxValue > checkMax

  end function any_vals_outside_2D
! ----------------------------------------------------------
  logical function any_vals_outside_3D(array, checkMin, checkMax)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: checkMin, checkMax

      ! Compact version using intrinsics below
      ! but an explicit loop is the only current solution on GPUs
    real(wp) :: minValue, maxValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    maxValue = maxval(array)
    !$acc end kernels
    any_vals_outside_3D = minValue < checkMin .or. maxValue > checkMax

  end function any_vals_outside_3D
  ! ----------------------------------------------------------
  ! Masked versions
  ! ----------------------------------------------------------
  logical function any_vals_outside_1D_masked(array, mask, checkMin, checkMax)
    real   (wp), dimension(:), intent(in) :: array
    logical(wl), dimension(:), intent(in) :: mask
    real(wp),                  intent(in) :: checkMin, checkMax

    real(wp) :: minValue, maxValue

    !$acc kernels copyin(array)
    minValue = minval(array, mask=mask)
    maxValue = maxval(array, mask=mask)
    !$acc end kernels
    any_vals_outside_1D_masked = minValue < checkMin .or. maxValue > checkMax

  end function any_vals_outside_1D_masked
! ----------------------------------------------------------
  logical function any_vals_outside_2D_masked(array, mask, checkMin, checkMax)
    real   (wp), dimension(:,:), intent(in) :: array
    logical(wl), dimension(:,:), intent(in) :: mask
    real(wp),                    intent(in) :: checkMin, checkMax

    real(wp) :: minValue, maxValue

    !$acc kernels copyin(array)
    minValue = minval(array, mask=mask)
    maxValue = maxval(array, mask=mask)
    !$acc end kernels
    any_vals_outside_2D_masked = minValue < checkMin .or. maxValue > checkMax

  end function any_vals_outside_2D_masked
! ----------------------------------------------------------
  logical function any_vals_outside_3D_masked(array, mask, checkMin, checkMax)
    real   (wp), dimension(:,:,:), intent(in) :: array
    logical(wl), dimension(:,:,:), intent(in) :: mask
    real(wp),                      intent(in) :: checkMin, checkMax

    real(wp) :: minValue, maxValue

    !$acc kernels copyin(array)
    minValue = minval(array, mask=mask)
    maxValue = maxval(array, mask=mask)
    !$acc end kernels
    any_vals_outside_3D_masked = minValue < checkMin .or. maxValue > checkMax

  end function any_vals_outside_3D_masked
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Extents
  !
  ! --------------------------------------------------------------------------------------
  function extents_are_1d(array, n1)
    real(wp), dimension(:          ), intent(in) :: array
    integer,                          intent(in) :: n1
    logical(wl)                                  :: extents_are_1d

    extents_are_1d = (size(array,1) == n1)
  end function extents_are_1d
  ! --------------------------------------------------------------------------------------
  function extents_are_2d(array, n1, n2)
    real(wp), dimension(:,:        ), intent(in) :: array
    integer,                          intent(in) :: n1, n2
    logical(wl)                                  :: extents_are_2d

    extents_are_2d = (size(array,1) == n1 .and. &
                      size(array,2) == n2 )
  end function extents_are_2d
  ! --------------------------------------------------------------------------------------
  function extents_are_3d(array, n1, n2, n3)
    real(wp), dimension(:,:,:      ), intent(in) :: array
    integer,                          intent(in) :: n1, n2, n3
    logical(wl)                                  :: extents_are_3d

    extents_are_3d = (size(array,1) == n1 .and. &
                      size(array,2) == n2 .and. &
                      size(array,3) == n3)
  end function extents_are_3d
  ! --------------------------------------------------------------------------------------
  function extents_are_4d(array, n1, n2, n3, n4)
    real(wp), dimension(:,:,:,:    ), intent(in) :: array
    integer,                          intent(in) :: n1, n2, n3, n4
    logical(wl)                                  :: extents_are_4d

    extents_are_4d = (size(array,1) == n1 .and. &
                      size(array,2) == n2 .and. &
                      size(array,3) == n3 .and. &
                      size(array,4) == n4)
  end function extents_are_4d
  ! --------------------------------------------------------------------------------------
  function extents_are_5d(array, n1, n2, n3, n4, n5)
    real(wp), dimension(:,:,:,:,:  ), intent(in) :: array
    integer,                          intent(in) :: n1, n2, n3, n4, n5
    logical(wl)                                  :: extents_are_5d

    extents_are_5d = (size(array,1) == n1 .and. &
                      size(array,2) == n2 .and. &
                      size(array,3) == n3 .and. &
                      size(array,4) == n4 .and. &
                      size(array,5) == n5 )
  end function extents_are_5d
  ! --------------------------------------------------------------------------------------
  function extents_are_6d(array, n1, n2, n3, n4, n5, n6)
    real(wp), dimension(:,:,:,:,:,:), intent(in) :: array
    integer,                          intent(in) :: n1, n2, n3, n4, n5, n6
    logical(wl)                                  :: extents_are_6d

    extents_are_6d = (size(array,1) == n1 .and. &
                      size(array,2) == n2 .and. &
                      size(array,3) == n3 .and. &
                      size(array,4) == n4 .and. &
                      size(array,5) == n5 .and. &
                      size(array,6) == n6 )
  end function extents_are_6d
  ! --------------------------------------------------------------------------------------
  function extents_are_2d_int(array, n1, n2)
    integer,  dimension(:,:        ), intent(in) :: array
    integer,                          intent(in) :: n1, n2
    logical(wl)                                  :: extents_are_2d_int

    extents_are_2d_int = (size(array,1) == n1 .and. &
                          size(array,2) == n2 )
  end function extents_are_2d_int
  !-------------------------------------------------------------------------------------------------
  ! Initializing arrays to 0
  !-------------------------------------------------------------------------------------------------
  subroutine zero_array_1D(ni, array) bind(C, name="zero_array_1D")
    integer,                 intent(in ) :: ni
    real(wp), dimension(ni), intent(out) :: array
    ! -----------------------
    integer :: i
    ! -----------------------
    !$acc parallel loop copyout(array)
    do i = 1, ni
      array(i) = 0.0_wp
    end do
  end subroutine zero_array_1D
  ! ----------------------------------------------------------
  subroutine zero_array_3D(ni, nj, nk, array) bind(C, name="zero_array_3D")
    integer,                         intent(in ) :: ni, nj, nk
    real(wp), dimension(ni, nj, nk), intent(out) :: array
    ! -----------------------
    integer :: i,j,k
    ! -----------------------
    !$acc parallel loop collapse(3) copyout(array)
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          array(i,j,k) = 0.0_wp
        end do
      end do
    end do

  end subroutine zero_array_3D
  ! ----------------------------------------------------------
  subroutine zero_array_4D(ni, nj, nk, nl, array) bind(C, name="zero_array_4D")
    integer,                             intent(in ) :: ni, nj, nk, nl
    real(wp), dimension(ni, nj, nk, nl), intent(out) :: array
    ! -----------------------
    integer :: i,j,k,l
    ! -----------------------
    !$acc parallel loop collapse(4) copyout(array)
    do l = 1, nl
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            array(i,j,k,l) = 0.0_wp
          end do
        end do
      end do
    end do

  end subroutine zero_array_4D

end module mo_rte_util_array
