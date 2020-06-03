! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Kernels to permute arrays

module mo_rrtmgp_util_reorder_kernels
  use mo_rte_kind,      only: wp
  implicit none
  public
contains
  ! ----------------------------------------------------------------------------
  subroutine reorder_123x312_kernel(d1, d2, d3, array_in, array_out) &
      bind(C, name = "reorder_123x312_kernel")
    integer,                         intent( in) :: d1, d2, d3
    real(wp), dimension(d1, d2, d3), intent( in) :: array_in
    real(wp), dimension(d3, d1, d2), intent(out) :: array_out

    integer :: i1, i2, i3, i10, i30, i1diff, i3diff
    integer, parameter :: tile = 32

    ! This kernel uses blocking to speed-up the transposition
    ! We read the data block by block (three outer loops)
    !  such that a block fits into fastest cache and the memory reads
    !  are resolved in the cache. The writes are contiguous here, so
    !  shouldn't be a problem.
    !  Tile size of 32x32 is empirical: big enough to read from the whole
    !  cache line, and small enough to fit into cache. Other numbers
    !  may give slightly better performance on different hardware.
    !
    !$acc parallel vector_length(tile*tile) &
    !$acc&     copyout(array_out) &
    !$acc&     copyin(array_in)
    !$acc loop gang collapse(3)
    do i2 = 1, d2
      do i10 = 1, d1, tile
        do i30 = 1, d3, tile

          !$acc loop vector collapse(2)
          do i1diff = 0, tile-1
            do i3diff = 0, tile-1
              i1 = i10 + i1diff
              i3 = i30 + i3diff
              if (i1 > d1 .or. i3 > d3) cycle

              array_out(i3,i1,i2) = array_in(i1,i2,i3)
            end do
          end do

        end do
      end do
    end do
    !$acc end parallel

  end subroutine reorder_123x312_kernel
  ! ----------------------------------------------------------------------------
  subroutine reorder_123x321_kernel(d1, d2, d3, array_in, array_out) &
      bind(C, name="reorder_123x321_kernel")
    integer,                         intent( in) :: d1, d2, d3
    real(wp), dimension(d1, d2, d3), intent( in) :: array_in
    real(wp), dimension(d3, d2, d1), intent(out) :: array_out

    integer :: i1, i2, i3, i10, i30, i1diff, i3diff, idiff
    integer, parameter :: tile = 32

    ! See the comment above
    !
    !$acc parallel vector_length(tile*tile) &
    !$acc&     copyout(array_out) &
    !$acc&     copyin(array_in)
    !$acc loop gang collapse(3)
    ! private(cache(:,:))
    do i2 = 1, d2
      do i10 = 1, d1, tile
        do i30 = 1, d3, tile

          !$acc loop vector collapse(2)
          do i1diff = 0, tile-1
            do i3diff = 0, tile-1
              i1 = i10 + i1diff
              i3 = i30 + i3diff
              if (i1 > d1 .or. i3 > d3) cycle

              array_out(i3,i2,i1) = array_in(i1,i2,i3)
            end do
          end do

        end do
      end do
    end do
    !$acc end parallel

  end subroutine reorder_123x321_kernel
  ! ----------------------------------------------------------------------------
end module mo_rrtmgp_util_reorder_kernels
