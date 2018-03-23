! This wrapper has been imported from PyCLES
! Add appropriate copyrights

module rrtmg_lw_wrapper

use iso_c_binding, only: c_double, c_int
use parrrtm, only : nbndlw
use rrtmg_lw_init, only: rrtmg_lw_ini
use rrtmg_lw_rad,  only: rrtmg_lw

implicit none

contains

subroutine c_rrtmg_lw_init(cpdair) bind(c)
    real(c_double), intent(in) :: cpdair
    call rrtmg_lw_ini(cpdair)
end subroutine c_rrtmg_lw_init


subroutine c_rrtmg_lw &
            (ncol    ,nlay    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
             inflglw ,iceflglw,liqflglw,cldfr   , &
             taucld  ,cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
             duflx_dt,duflxc_dt ) bind(c)
      integer(c_int), intent(in) :: ncol            ! Number of horizontal columns
      integer(c_int), intent(in) :: nlay            ! Number of model layers
      integer(c_int), intent(inout) :: icld         ! Cloud overlap method
                                                      !    0: Clear only
                                                      !    1: Random
                                                      !    2: Maximum/random
                                                      !    3: Maximum
      integer(c_int), intent(in) :: idrv            ! Flag for calculation of dFdT, the change
                                                      !    in upward flux as a function of 
                                                      !    surface temperature [0=off, 1=on]
                                                      !    0: Normal forward calculation
                                                      !    1: Normal forward calculation with
                                                      !       duflx_dt and duflxc_dt output

      real(c_double), intent(in) :: play(ncol,nlay)          ! Layer pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: plev(ncol,nlay+1)          ! Interface pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay+1)
      real(c_double), intent(in) :: tlay(ncol,nlay)          ! Layer temperatures (K)
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: tlev(ncol,nlay+1)          ! Interface temperatures (K)
                                                      !    Dimensions: (ncol,nlay+1)
      real(c_double), intent(in) :: tsfc(ncol)            ! Surface temperature (K)
                                                      !    Dimensions: (ncol)
      real(c_double), intent(in) :: h2ovmr(ncol,nlay)        ! H2O volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: o3vmr(ncol,nlay)         ! O3 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: co2vmr(ncol,nlay)        ! CO2 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: ch4vmr(ncol,nlay)        ! Methane volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: n2ovmr(ncol,nlay)        ! Nitrous oxide volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: o2vmr(ncol,nlay)         ! Oxygen volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: cfc11vmr(ncol,nlay)      ! CFC11 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: cfc12vmr(ncol,nlay)      ! CFC12 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: cfc22vmr(ncol,nlay)      ! CFC22 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: ccl4vmr(ncol,nlay)       ! CCL4 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: emis(ncol,nbndlw)          ! Surface emissivity
                                                      !    Dimensions: (ncol,nbndlw)

      integer(c_int), intent(in) :: inflglw         ! Flag for cloud optical properties
      integer(c_int), intent(in) :: iceflglw        ! Flag for ice particle specification
      integer(c_int), intent(in) :: liqflglw        ! Flag for liquid droplet specification

      real(c_double), intent(in) :: cldfr(ncol,nlay)         ! Cloud fraction
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: cicewp(ncol,nlay)        ! Cloud ice water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: cliqwp(ncol,nlay)        ! Cloud liquid water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: reice(ncol,nlay)         ! Cloud ice particle effective size (microns)
                                                      !    Dimensions: (ncol,nlay)
                                                      ! specific definition of reice depends on setting of iceflglw:
                                                      ! iceflglw = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !               r_ec must be >= 10.0 microns
                                                      ! iceflglw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !               r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflglw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !               r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflglw = 3: generalized effective size, dge, (Fu, 1996),
                                                      !               dge range is limited to 5.0 to 140.0 microns
                                                      !               [dge = 1.0315 * r_ec]
      real(c_double), intent(in) :: reliq(ncol,nlay)         ! Cloud water drop effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(in) :: taucld(nbndlw,ncol,nlay)      ! In-cloud optical depth
                                                      !    Dimensions: (nbndlw,ncol,nlay)
      real(c_double), intent(in) :: tauaer(ncol,nlay,nbndlw)      ! aerosol optical depth
                                                      !    Dimensions: (ncol,nlay,nbndlw)


! ----- Output -----

      real(c_double), intent(out) :: uflx(ncol,nlay+1)         ! Total sky longwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(c_double), intent(out) :: dflx(ncol,nlay+1)         ! Total sky longwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(c_double), intent(out) :: hr(ncol,nlay)           ! Total sky longwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)
      real(c_double), intent(out) :: uflxc(ncol,nlay+1)        ! Clear sky longwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(c_double), intent(out) :: dflxc(ncol,nlay+1)        ! Clear sky longwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(c_double), intent(out) :: hrc(ncol,nlay)          ! Clear sky longwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)

! ----- Optional Output -----
      real(c_double), intent(out) :: duflx_dt(ncol,nlay+1)     
                                                      ! change in upward longwave flux (w/m2/k)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (ncol,nlay+1)
      real(c_double), intent(out) :: duflxc_dt(ncol,nlay+1)    
                                                      ! change in clear sky upward longwave flux (w/m2/k)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (ncol,nlay+1)
    
    
    call rrtmg_lw &
            (ncol    ,nlay    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
             inflglw ,iceflglw,liqflglw,cldfr   , &
             taucld  ,cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
             duflx_dt,duflxc_dt )
end subroutine c_rrtmg_lw


end module
