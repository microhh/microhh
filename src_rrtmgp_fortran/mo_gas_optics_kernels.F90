! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Numeric calculations for gas optics. Absorption and Rayleigh optical depths,
!   source functions.

module mo_gas_optics_kernels
  use mo_rte_kind,      only : wp, wl
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  ! Compute interpolation coefficients
  ! for calculations of major optical depths, minor optical depths, Rayleigh,
  ! and Planck fractions
  subroutine interpolation( &
                ncol,nlay,ngas,nflav,neta, npres, ntemp, &
                flavor,                                  &
                press_ref_log, temp_ref,press_ref_log_delta,    &
                temp_ref_min,temp_ref_delta,press_ref_trop_log, &
                vmr_ref,                                        &
                play,tlay,col_gas,                              &
                jtemp,fmajor,fminor,col_mix,tropo,jeta,jpress) bind(C, name="interpolation")
    ! input dimensions
    integer,                            intent(in) :: ncol,nlay
    integer,                            intent(in) :: ngas,nflav,neta,npres,ntemp
    integer,     dimension(2,nflav),    intent(in) :: flavor
    real(wp),    dimension(npres),      intent(in) :: press_ref_log
    real(wp),    dimension(ntemp),      intent(in) :: temp_ref
    real(wp),                           intent(in) :: press_ref_log_delta, &
                                                      temp_ref_min, temp_ref_delta, &
                                                      press_ref_trop_log
    real(wp),    dimension(2,0:ngas,ntemp), intent(in) :: vmr_ref

    ! inputs from profile or parent function
    real(wp),    dimension(ncol,nlay),        intent(in) :: play, tlay
    real(wp),    dimension(ncol,nlay,0:ngas), intent(in) :: col_gas

    ! outputs
    integer,     dimension(ncol,nlay), intent(out) :: jtemp, jpress
    logical(wl), dimension(ncol,nlay), intent(out) :: tropo
    integer,     dimension(2,    nflav,ncol,nlay), intent(out) :: jeta
    real(wp),    dimension(2,    nflav,ncol,nlay), intent(out) :: col_mix
    real(wp),    dimension(2,2,2,nflav,ncol,nlay), intent(out) :: fmajor
    real(wp),    dimension(2,2,  nflav,ncol,nlay), intent(out) :: fminor
    ! -----------------
    ! local
    real(wp), dimension(ncol,nlay) :: ftemp, fpress ! interpolation fraction for temperature, pressure
    real(wp) :: locpress ! needed to find location in pressure grid
    real(wp) :: ratio_eta_half ! ratio of vmrs of major species that defines eta=0.5
                               ! for given flavor and reference temperature level
    real(wp) :: eta, feta      ! binary_species_parameter, interpolation variable for eta
    real(wp) :: loceta         ! needed to find location in eta grid
    real(wp) :: ftemp_term
    ! -----------------
    ! local indexes
    integer :: icol, ilay, iflav, igases(2), itropo, itemp

    do ilay = 1, nlay
      do icol = 1, ncol
        ! index and factor for temperature interpolation
        jtemp(icol,ilay) = int((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta)
        jtemp(icol,ilay) = min(ntemp - 1, max(1, jtemp(icol,ilay))) ! limit the index range
        ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta

        ! index and factor for pressure interpolation
        locpress = 1._wp + (log(play(icol,ilay)) - press_ref_log(1)) / press_ref_log_delta
        jpress(icol,ilay) = min(npres-1, max(1, int(locpress)))
        fpress(icol,ilay) = locpress - float(jpress(icol,ilay))

        ! determine if in lower or upper part of atmosphere
        tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log
      end do
    end do

    do ilay = 1, nlay
      do icol = 1, ncol
        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))
        ! loop over implemented combinations of major species
        do iflav = 1, nflav
          igases(:) = flavor(:,iflav)
          do itemp = 1, 2
            ! compute interpolation fractions needed for lower, then upper reference temperature level
            ! compute binary species parameter (eta) for flavor and temperature and
            !  associated interpolation index and factors
            ratio_eta_half = vmr_ref(itropo,igases(1),(jtemp(icol,ilay)+itemp-1)) / &
                             vmr_ref(itropo,igases(2),(jtemp(icol,ilay)+itemp-1))
            col_mix(itemp,iflav,icol,ilay) = col_gas(icol,ilay,igases(1)) + ratio_eta_half * col_gas(icol,ilay,igases(2))
            eta = merge(col_gas(icol,ilay,igases(1)) / col_mix(itemp,iflav,icol,ilay), 0.5_wp, &
                        col_mix(itemp,iflav,icol,ilay) > 2._wp * tiny(col_mix))
            loceta = eta * float(neta-1)
            jeta(itemp,iflav,icol,ilay) = min(int(loceta)+1, neta-1)
            feta = mod(loceta, 1.0_wp)
            ! compute interpolation fractions needed for minor species
            ! ftemp_term = (1._wp-ftemp(icol,ilay)) for itemp = 1, ftemp(icol,ilay) for itemp=2
            ftemp_term = (real(2-itemp, wp) + real(2*itemp-3, wp) * ftemp(icol,ilay))
            fminor(1,itemp,iflav,icol,ilay) = (1._wp-feta) * ftemp_term
            fminor(2,itemp,iflav,icol,ilay) =        feta  * ftemp_term
            ! compute interpolation fractions needed for major species
            fmajor(1,1,itemp,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(1,itemp,iflav,icol,ilay)
            fmajor(2,1,itemp,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(2,itemp,iflav,icol,ilay)
            fmajor(1,2,itemp,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(1,itemp,iflav,icol,ilay)
            fmajor(2,2,itemp,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(2,itemp,iflav,icol,ilay)
          end do ! reference temperatures
        end do ! iflav
      end do ! icol,ilay
    end do

  end subroutine interpolation
  ! --------------------------------------------------------------------------------------
  !
  ! Compute minor and major species opitcal depth from pre-computed interpolation coefficients
  !   (jeta,jtemp,jpress)
  !
  subroutine compute_tau_absorption(                &
                ncol,nlay,nbnd,ngpt,                &  ! dimensions
                ngas,nflav,neta,npres,ntemp,        &
                nminorlower, nminorklower,          & ! number of minor contributors, total num absorption coeffs
                nminorupper, nminorkupper,          &
                idx_h2o,                            &
                gpoint_flavor,                      &
                band_lims_gpt,                      &
                kmajor,                             &
                kminor_lower,                       &
                kminor_upper,                       &
                minor_limits_gpt_lower,             &
                minor_limits_gpt_upper,             &
                minor_scales_with_density_lower,    &
                minor_scales_with_density_upper,    &
                scale_by_complement_lower,          &
                scale_by_complement_upper,          &
                idx_minor_lower,                    &
                idx_minor_upper,                    &
                idx_minor_scaling_lower,            &
                idx_minor_scaling_upper,            &
                kminor_start_lower,                 &
                kminor_start_upper,                 &
                tropo,                              &
                col_mix,fmajor,fminor,              &
                play,tlay,col_gas,                  &
                jeta,jtemp,jpress,                  &
                tau) bind(C, name="compute_tau_absorption")
    ! ---------------------
    ! input dimensions
    integer,                                intent(in) :: ncol,nlay,nbnd,ngpt
    integer,                                intent(in) :: ngas,nflav,neta,npres,ntemp
    integer,                                intent(in) :: nminorlower, nminorklower,nminorupper, nminorkupper
    integer,                                intent(in) :: idx_h2o
    ! ---------------------
    ! inputs from object
    integer,     dimension(2,ngpt),                  intent(in) :: gpoint_flavor
    integer,     dimension(2,nbnd),                  intent(in) :: band_lims_gpt
    real(wp),    dimension(ngpt,neta,npres+1,ntemp), intent(in) :: kmajor
    real(wp),    dimension(nminorklower,neta,ntemp), intent(in) :: kminor_lower
    real(wp),    dimension(nminorkupper,neta,ntemp), intent(in) :: kminor_upper
    integer,     dimension(2,nminorlower),           intent(in) :: minor_limits_gpt_lower
    integer,     dimension(2,nminorupper),           intent(in) :: minor_limits_gpt_upper
    logical(wl), dimension(  nminorlower),           intent(in) :: minor_scales_with_density_lower
    logical(wl), dimension(  nminorupper),           intent(in) :: minor_scales_with_density_upper
    logical(wl), dimension(  nminorlower),           intent(in) :: scale_by_complement_lower
    logical(wl), dimension(  nminorupper),           intent(in) :: scale_by_complement_upper
    integer,     dimension(  nminorlower),           intent(in) :: idx_minor_lower
    integer,     dimension(  nminorupper),           intent(in) :: idx_minor_upper
    integer,     dimension(  nminorlower),           intent(in) :: idx_minor_scaling_lower
    integer,     dimension(  nminorupper),           intent(in) :: idx_minor_scaling_upper
    integer,     dimension(  nminorlower),           intent(in) :: kminor_start_lower
    integer,     dimension(  nminorupper),           intent(in) :: kminor_start_upper
    logical(wl), dimension(ncol,nlay),               intent(in) :: tropo
    ! ---------------------
    ! inputs from profile or parent function
    real(wp), dimension(2,    nflav,ncol,nlay       ), intent(in) :: col_mix
    real(wp), dimension(2,2,2,nflav,ncol,nlay       ), intent(in) :: fmajor
    real(wp), dimension(2,2,  nflav,ncol,nlay       ), intent(in) :: fminor
    real(wp), dimension(            ncol,nlay       ), intent(in) :: play, tlay      ! pressure and temperature
    real(wp), dimension(            ncol,nlay,0:ngas), intent(in) :: col_gas
    integer,  dimension(2,    nflav,ncol,nlay       ), intent(in) :: jeta
    integer,  dimension(            ncol,nlay       ), intent(in) :: jtemp
    integer,  dimension(            ncol,nlay       ), intent(in) :: jpress
    ! ---------------------
    ! output - optical depth
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! ---------------------
    ! Local variables
    !
    logical                    :: top_at_1
    integer, dimension(ncol,2) :: itropo_lower, itropo_upper
    ! ----------------------------------------------------------------

    ! ---------------------
    ! Layer limits of upper, lower atmospheres
    ! ---------------------
    top_at_1 = play(1,1) < play(1, nlay)
    if(top_at_1) then
      itropo_lower(:, 1) = minloc(play, dim=2, mask=tropo)
      itropo_lower(:, 2) = nlay
      itropo_upper(:, 1) = 1
      itropo_upper(:, 2) = maxloc(play, dim=2, mask=(.not. tropo))
    else
      itropo_lower(:, 1) = 1
      itropo_lower(:, 2) = minloc(play, dim=2, mask= tropo)
      itropo_upper(:, 1) = maxloc(play, dim=2, mask=(.not. tropo))
      itropo_upper(:, 2) = nlay
    end if
    ! ---------------------
    ! Major Species
    ! ---------------------
    call gas_optical_depths_major(   &
          ncol,nlay,nbnd,ngpt,       & ! dimensions
          nflav,neta,npres,ntemp,    &
          gpoint_flavor,             &
          band_lims_gpt,             &
          kmajor,                    &
          col_mix,fmajor,            &
          jeta,tropo,jtemp,jpress,   &
          tau)
    ! ---------------------
    ! Minor Species - lower
    ! ---------------------
    call gas_optical_depths_minor(     &
           ncol,nlay,ngpt,             & ! dimensions
           ngas,nflav,ntemp,neta,      &
           nminorlower,nminorklower,   &
           idx_h2o,                    &
           gpoint_flavor(1,:),         &
           kminor_lower,               &
           minor_limits_gpt_lower,     &
           minor_scales_with_density_lower, &
           scale_by_complement_lower,  &
           idx_minor_lower,            &
           idx_minor_scaling_lower,    &
           kminor_start_lower,         &
           play, tlay,                 &
           col_gas,fminor,jeta,        &
           itropo_lower,jtemp,         &
           tau)
    ! ---------------------
    ! Minor Species - upper
    ! ---------------------
    call gas_optical_depths_minor(     &
           ncol,nlay,ngpt,             & ! dimensions
           ngas,nflav,ntemp,neta,      &
           nminorupper,nminorkupper,   &
           idx_h2o,                    &
           gpoint_flavor(2,:),         &
           kminor_upper,               &
           minor_limits_gpt_upper,     &
           minor_scales_with_density_upper, &
           scale_by_complement_upper,  &
           idx_minor_upper,            &
           idx_minor_scaling_upper,    &
           kminor_start_upper,         &
           play, tlay,                 &
           col_gas,fminor,jeta,        &
           itropo_upper,jtemp,         &
           tau)
  end subroutine compute_tau_absorption
  ! --------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------
  !
  ! compute minor species optical depths
  !
  subroutine gas_optical_depths_major(ncol,nlay,nbnd,ngpt,&
                                      nflav,neta,npres,ntemp,      & ! dimensions
                                      gpoint_flavor, band_lims_gpt,   & ! inputs from object
                                      kmajor,                         &
                                      col_mix,fmajor,                 &
                                      jeta,tropo,jtemp,jpress,        & ! local input
                                      tau) bind(C, name="gas_optical_depths_major")
    ! input dimensions
    integer, intent(in) :: ncol, nlay, nbnd, ngpt, nflav,neta,npres,ntemp  ! dimensions

    ! inputs from object
    integer,  dimension(2,ngpt),  intent(in) :: gpoint_flavor
    integer,  dimension(2,nbnd),  intent(in) :: band_lims_gpt ! start and end g-point for each band
    real(wp), dimension(ngpt,neta,npres+1,ntemp), intent(in) :: kmajor

    ! inputs from profile or parent function
    real(wp),    dimension(2,    nflav,ncol,nlay), intent(in) :: col_mix
    real(wp),    dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
    integer,     dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
    logical(wl), dimension(ncol,nlay), intent(in) :: tropo
    integer,     dimension(ncol,nlay), intent(in) :: jtemp, jpress

    ! outputs
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! -----------------
    ! local variables
    real(wp) :: tau_major(ngpt) ! major species optical depth
    ! local index
    integer :: icol, ilay, iflav, ibnd, igpt, itropo
    integer :: gptS, gptE

    ! -----------------

    do icol = 1, ncol
      do ilay = 1, nlay
        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))
        ! optical depth calculation for major species
        do ibnd = 1, nbnd
          gptS = band_lims_gpt(1, ibnd)
          gptE = band_lims_gpt(2, ibnd)
          iflav = gpoint_flavor(itropo, gptS) !eta interpolation depends on band's flavor
          tau_major(gptS:gptE) = &
            ! interpolation in temperature, pressure, and eta
            interpolate3D_byflav(col_mix(:,iflav,icol,ilay),                                     &
                                 fmajor(:,:,:,iflav,icol,ilay), kmajor,                          &
                                 band_lims_gpt(1, ibnd), band_lims_gpt(2, ibnd),                 &
                                 jeta(:,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo)
          tau(gptS:gptE,ilay,icol) = tau(gptS:gptE,ilay,icol) + tau_major(gptS:gptE)
        end do ! igpt
      end do
    end do ! ilay
  end subroutine gas_optical_depths_major

  ! ----------------------------------------------------------
  !
  ! compute minor species optical depths
  !
  subroutine gas_optical_depths_minor(ncol,nlay,ngpt,        &
                                      ngas,nflav,ntemp,neta, &
                                      nminor,nminork,        &
                                      idx_h2o,               &
                                      gpt_flv,               &
                                      kminor,                &
                                      minor_limits_gpt,      &
                                      minor_scales_with_density,    &
                                      scale_by_complement,   &
                                      idx_minor, idx_minor_scaling, &
                                      kminor_start,          &
                                      play, tlay,            &
                                      col_gas,fminor,jeta,   &
                                      layer_limits,jtemp,    &
                                      tau) bind(C, name="gas_optical_depths_minor")
    integer,                                     intent(in   ) :: ncol,nlay,ngpt
    integer,                                     intent(in   ) :: ngas,nflav
    integer,                                     intent(in   ) :: ntemp,neta,nminor,nminork
    integer,                                     intent(in   ) :: idx_h2o
    integer,     dimension(ngpt),                intent(in   ) :: gpt_flv
    real(wp),    dimension(nminork,neta,ntemp),  intent(in   ) :: kminor
    integer,     dimension(2,nminor),            intent(in   ) :: minor_limits_gpt
    logical(wl), dimension(  nminor),            intent(in   ) :: minor_scales_with_density
    logical(wl), dimension(  nminor),            intent(in   ) :: scale_by_complement
    integer,     dimension(  nminor),            intent(in   ) :: kminor_start
    integer,     dimension(  nminor),            intent(in   ) :: idx_minor, idx_minor_scaling
    real(wp),    dimension(ncol,nlay),           intent(in   ) :: play, tlay
    real(wp),    dimension(ncol,nlay,0:ngas),    intent(in   ) :: col_gas
    real(wp),    dimension(2,2,nflav,ncol,nlay), intent(in   ) :: fminor
    integer,     dimension(2,  nflav,ncol,nlay), intent(in   ) :: jeta
    integer,     dimension(ncol, 2),             intent(in   ) :: layer_limits
    integer,     dimension(ncol,nlay),           intent(in   ) :: jtemp
    real(wp),    dimension(ngpt,nlay,ncol),      intent(inout) :: tau
    ! -----------------
    ! local variables
    real(wp), parameter :: PaTohPa = 0.01
    real(wp) :: vmr_fact, dry_fact             ! conversion from column abundance to dry vol. mixing ratio;
    real(wp) :: scaling, kminor_loc            ! minor species absorption coefficient, optical depth
    integer  :: icol, ilay, iflav, igpt, imnr
    integer  :: gptS, gptE
    real(wp), dimension(ngpt) :: tau_minor
    ! -----------------
    !
    ! Guard against layer limits being 0 -- that means don't do anything i.e. there are no
    !   layers with pressures in the upper or lower atmosphere respectively
    ! First check skips the routine entirely if all columns are out of bounds...
    !
    if(any(layer_limits(:,1) > 0)) then
      do imnr = 1, size(scale_by_complement,dim=1) ! loop over minor absorbers in each band
        do icol = 1, ncol
          !
          ! This check skips individual columns with no pressures in range
          !
          if(layer_limits(icol,1) > 0) then
            do ilay = layer_limits(icol,1), layer_limits(icol,2)
              !
              ! Scaling of minor gas absortion coefficient begins with column amount of minor gas
              !
              scaling = col_gas(icol,ilay,idx_minor(imnr))
              !
              ! Density scaling (e.g. for h2o continuum, collision-induced absorption)
              !
              if (minor_scales_with_density(imnr)) then
                !
                ! NOTE: P needed in hPa to properly handle density scaling.
                !
                scaling = scaling * (PaTohPa*play(icol,ilay)/tlay(icol,ilay))
                if(idx_minor_scaling(imnr) > 0) then  ! there is a second gas that affects this gas's absorption
                  vmr_fact = 1._wp / col_gas(icol,ilay,0)
                  dry_fact = 1._wp / (1._wp + col_gas(icol,ilay,idx_h2o) * vmr_fact)
                  ! scale by density of special gas
                  if (scale_by_complement(imnr)) then ! scale by densities of all gases but the special one
                    scaling = scaling * (1._wp - col_gas(icol,ilay,idx_minor_scaling(imnr)) * vmr_fact * dry_fact)
                  else
                    scaling = scaling *          col_gas(icol,ilay,idx_minor_scaling(imnr)) * vmr_fact * dry_fact
                  endif
                endif
              endif
              !
              ! Interpolation of absorption coefficient and calculation of optical depth
              !
              ! Which gpoint range does this minor gas affect?
              gptS = minor_limits_gpt(1,imnr)
              gptE = minor_limits_gpt(2,imnr)
              iflav = gpt_flv(gptS)
              tau_minor(gptS:gptE) = scaling *                   &
                                      interpolate2D_byflav(fminor(:,:,iflav,icol,ilay), &
                                                           kminor, &
                                                           kminor_start(imnr), kminor_start(imnr)+(gptE-gptS), &
                                                           jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
              tau(gptS:gptE,ilay,icol) = tau(gptS:gptE,ilay,icol) + tau_minor(gptS:gptE)
            enddo
          end if
        enddo
      enddo
    end if
  end subroutine gas_optical_depths_minor
  ! ----------------------------------------------------------
  !
  ! compute Rayleigh scattering optical depths
  !
  subroutine compute_tau_rayleigh(ncol,nlay,nbnd,ngpt,         &
                                  ngas,nflav,neta,npres,ntemp, &
                                  gpoint_flavor,band_lims_gpt, &
                                  krayl,                       &
                                  idx_h2o, col_dry,col_gas,    &
                                  fminor,jeta,tropo,jtemp,     &
                                  tau_rayleigh) bind(C, name="compute_tau_rayleigh")
    integer,                                     intent(in ) :: ncol,nlay,nbnd,ngpt
    integer,                                     intent(in ) :: ngas,nflav,neta,npres,ntemp
    integer,     dimension(2,ngpt),              intent(in ) :: gpoint_flavor
    integer,     dimension(2,nbnd),              intent(in ) :: band_lims_gpt ! start and end g-point for each band
    real(wp),    dimension(ngpt,neta,ntemp,2),   intent(in ) :: krayl
    integer,                                     intent(in ) :: idx_h2o
    real(wp),    dimension(ncol,nlay),           intent(in ) :: col_dry
    real(wp),    dimension(ncol,nlay,0:ngas),    intent(in ) :: col_gas
    real(wp),    dimension(2,2,nflav,ncol,nlay), intent(in ) :: fminor
    integer,     dimension(2,  nflav,ncol,nlay), intent(in ) :: jeta
    logical(wl), dimension(ncol,nlay),           intent(in ) :: tropo
    integer,     dimension(ncol,nlay),           intent(in ) :: jtemp
    ! outputs
    real(wp),    dimension(ngpt,nlay,ncol),      intent(out) :: tau_rayleigh
    ! -----------------
    ! local variables
    real(wp) :: k(ngpt) ! rayleigh scattering coefficient
    integer  :: icol, ilay, iflav, ibnd, igpt, gptS, gptE
    integer  :: itropo
    ! -----------------
    do ilay = 1, nlay
      do icol = 1, ncol
        itropo = merge(1,2,tropo(icol,ilay)) ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        do ibnd = 1, nbnd
          gptS = band_lims_gpt(1, ibnd)
          gptE = band_lims_gpt(2, ibnd)
          iflav = gpoint_flavor(itropo, gptS) !eta interpolation depends on band's flavor
          k(gptS:gptE) = interpolate2D_byflav(fminor(:,:,iflav,icol,ilay), &
                                              krayl(:,:,:,itropo),      &
                                              gptS, gptE, jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
          tau_rayleigh(gptS:gptE,ilay,icol) = k(gptS:gptE) * &
                                              (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay))
        end do
      end do
    end do
  end subroutine compute_tau_rayleigh

  ! ----------------------------------------------------------
  subroutine compute_Planck_source(                        &
                    ncol, nlay, nbnd, ngpt,                &
                    nflav, neta, npres, ntemp, nPlanckTemp,&
                    tlay, tlev, tsfc, sfc_lay,             &
                    fmajor, jeta, tropo, jtemp, jpress,    &
                    gpoint_bands, band_lims_gpt,           &
                    pfracin, temp_ref_min, totplnk_delta, totplnk, gpoint_flavor, &
                    sfc_src, lay_src, lev_src_inc, lev_src_dec) bind(C, name="compute_Planck_source")
    integer,                                    intent(in) :: ncol, nlay, nbnd, ngpt
    integer,                                    intent(in) :: nflav, neta, npres, ntemp, nPlanckTemp
    real(wp),    dimension(ncol,nlay  ),        intent(in) :: tlay
    real(wp),    dimension(ncol,nlay+1),        intent(in) :: tlev
    real(wp),    dimension(ncol       ),        intent(in) :: tsfc
    integer,                                    intent(in) :: sfc_lay
    ! Interpolation variables
    real(wp),    dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
    integer,     dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
    logical(wl), dimension(            ncol,nlay), intent(in) :: tropo
    integer,     dimension(            ncol,nlay), intent(in) :: jtemp, jpress
    ! Table-specific
    integer, dimension(ngpt),                     intent(in) :: gpoint_bands ! start and end g-point for each band
    integer, dimension(2, nbnd),                  intent(in) :: band_lims_gpt ! start and end g-point for each band
    real(wp),                                     intent(in) :: temp_ref_min, totplnk_delta
    real(wp), dimension(ngpt,neta,npres+1,ntemp), intent(in) :: pfracin
    real(wp), dimension(nPlanckTemp,nbnd),        intent(in) :: totplnk
    integer,  dimension(2,ngpt),                  intent(in) :: gpoint_flavor

    real(wp), dimension(ngpt,     ncol), intent(out) :: sfc_src
    real(wp), dimension(ngpt,nlay,ncol), intent(out) :: lay_src
    real(wp), dimension(ngpt,nlay,ncol), intent(out) :: lev_src_inc, lev_src_dec
    ! -----------------
    ! local
    integer  :: ilay, icol, igpt, ibnd, itropo, iflav
    integer  :: gptS, gptE
    real(wp), dimension(2), parameter :: one = [1._wp, 1._wp]
    real(wp) :: pfrac          (ngpt,nlay,  ncol)
    real(wp) :: planck_function(nbnd,nlay+1,ncol)
    ! -----------------

    ! Calculation of fraction of band's Planck irradiance associated with each g-point
    do icol = 1, ncol
      do ilay = 1, nlay
        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))
        do ibnd = 1, nbnd
          gptS = band_lims_gpt(1, ibnd)
          gptE = band_lims_gpt(2, ibnd)
          iflav = gpoint_flavor(itropo, gptS) !eta interpolation depends on band's flavor
          pfrac(gptS:gptE,ilay,icol) = &
            ! interpolation in temperature, pressure, and eta
            interpolate3D_byflav(one, fmajor(:,:,:,iflav,icol,ilay), pfracin, &
                          band_lims_gpt(1, ibnd), band_lims_gpt(2, ibnd),                 &
                          jeta(:,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo)
        end do ! band
      end do   ! layer
    end do     ! column

    !
    ! Planck function by band for the surface
    ! Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
    !
    do icol = 1, ncol
      planck_function(1:nbnd,1,icol) = interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk)
      !
      ! Map to g-points
      !
      do ibnd = 1, nbnd
        gptS = band_lims_gpt(1, ibnd)
        gptE = band_lims_gpt(2, ibnd)
        do igpt = gptS, gptE
          sfc_src(igpt, icol) = pfrac(igpt,sfc_lay,icol) * planck_function(ibnd, 1, icol)
        end do
      end do
    end do ! icol

    do icol = 1, ncol
      do ilay = 1, nlay
        ! Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
        planck_function(1:nbnd,ilay,icol) = interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk)
        !
        ! Map to g-points
        !
        do ibnd = 1, nbnd
          gptS = band_lims_gpt(1, ibnd)
          gptE = band_lims_gpt(2, ibnd)
          do igpt = gptS, gptE
            lay_src(igpt,ilay,icol) = pfrac(igpt,ilay,icol) * planck_function(ibnd,ilay,icol)
          end do
        end do
      end do ! ilay
    end do ! icol

    ! compute level source irradiances for each g-point, one each for upward and downward paths
    do icol = 1, ncol
      planck_function(1:nbnd,       1,icol) = interpolate1D(tlev(icol,     1), temp_ref_min, totplnk_delta, totplnk)
      do ilay = 1, nlay
        planck_function(1:nbnd,ilay+1,icol) = interpolate1D(tlev(icol,ilay+1), temp_ref_min, totplnk_delta, totplnk)
        !
        ! Map to g-points
        !
        do ibnd = 1, nbnd
          gptS = band_lims_gpt(1, ibnd)
          gptE = band_lims_gpt(2, ibnd)
          do igpt = gptS, gptE
            lev_src_inc(igpt,ilay,icol) = pfrac(igpt,ilay,icol) * planck_function(ibnd,ilay+1,icol)
            lev_src_dec(igpt,ilay,icol) = pfrac(igpt,ilay,icol) * planck_function(ibnd,ilay,  icol)
          end do
        end do
      end do ! ilay
    end do ! icol

  end subroutine compute_Planck_source
  ! ----------------------------------------------------------
  !
  ! One dimensional interpolation -- return all values along second table dimension
  !
  pure function interpolate1D(val, offset, delta, table) result(res)
    ! input
    real(wp), intent(in) :: val,    & ! axis value at which to evaluate table
                            offset, & ! minimum of table axis
                            delta     ! step size of table axis
    real(wp), dimension(:,:), &
              intent(in) :: table ! dimensions (axis, values)
    ! output
    real(wp), dimension(size(table,dim=2)) :: res

    ! local
    real(wp) :: val0 ! fraction index adjusted by offset and delta
    integer :: index ! index term
    real(wp) :: frac ! fractional term
    ! -------------------------------------
    val0 = (val - offset) / delta
    frac = val0 - int(val0) ! get fractional part
    index = min(size(table,dim=1)-1, max(1, int(val0)+1)) ! limit the index range
    res(:) = table(index,:) + frac * (table(index+1,:) - table(index,:))
  end function interpolate1D
  ! ----------------------------------------------------------------------------------------
  !   This function returns a single value from a subset (in gpoint) of the k table
  !
  pure function interpolate2D(fminor, k, igpt, jeta, jtemp) result(res)
    real(wp), dimension(2,2), intent(in) :: fminor ! interpolation fractions for minor species
                                       ! index(1) : reference eta level (temperature dependent)
                                       ! index(2) : reference temperature level
    real(wp), dimension(:,:,:), intent(in) :: k ! (g-point, eta, temp)
    integer,                    intent(in) :: igpt, jtemp ! interpolation index for temperature
    integer, dimension(2),      intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    real(wp)                             :: res ! the result

    res =  &
      fminor(1,1) * k(igpt, jeta(1)  , jtemp  ) + &
      fminor(2,1) * k(igpt, jeta(1)+1, jtemp  ) + &
      fminor(1,2) * k(igpt, jeta(2)  , jtemp+1) + &
      fminor(2,2) * k(igpt, jeta(2)+1, jtemp+1)
  end function interpolate2D
  ! ----------------------------------------------------------
  !   This function returns a range of values from a subset (in gpoint) of the k table
  !
  pure function interpolate2D_byflav(fminor, k, gptS, gptE, jeta, jtemp) result(res)
    real(wp), dimension(2,2), intent(in) :: fminor ! interpolation fractions for minor species
                                       ! index(1) : reference eta level (temperature dependent)
                                       ! index(2) : reference temperature level
    real(wp), dimension(:,:,:), intent(in) :: k ! (g-point, eta, temp)
    integer,                    intent(in) :: gptS, gptE, jtemp ! interpolation index for temperature
    integer, dimension(2),      intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    real(wp), dimension(gptE-gptS+1)       :: res ! the result

    ! Local variable
    integer :: igpt
    ! each code block is for a different reference temperature
    do igpt = 1, gptE-gptS+1
      res(igpt) = fminor(1,1) * k(gptS+igpt-1, jeta(1)  , jtemp  ) + &
                  fminor(2,1) * k(gptS+igpt-1, jeta(1)+1, jtemp  ) + &
                  fminor(1,2) * k(gptS+igpt-1, jeta(2)  , jtemp+1) + &
                  fminor(2,2) * k(gptS+igpt-1, jeta(2)+1, jtemp+1)
    end do
  end function interpolate2D_byflav
  ! ----------------------------------------------------------
  ! interpolation in temperature, pressure, and eta
  pure function interpolate3D(scaling, fmajor, k, igpt, jeta, jtemp, jpress) result(res)
    real(wp), dimension(2),     intent(in) :: scaling
    real(wp), dimension(2,2,2), intent(in) :: fmajor ! interpolation fractions for major species
                                                     ! index(1) : reference eta level (temperature dependent)
                                                     ! index(2) : reference pressure level
                                                     ! index(3) : reference temperature level
    real(wp), dimension(:,:,:,:),intent(in) :: k ! (gpt, eta,temp,press)
    integer,                     intent(in) :: igpt
    integer, dimension(2),       intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    integer,                     intent(in) :: jtemp ! interpolation index for temperature
    integer,                     intent(in) :: jpress ! interpolation index for pressure
    real(wp)                                :: res ! the result
    ! each code block is for a different reference temperature
    res =  &
      scaling(1) * &
      ( fmajor(1,1,1) * k(igpt, jeta(1)  , jpress-1, jtemp  ) + &
        fmajor(2,1,1) * k(igpt, jeta(1)+1, jpress-1, jtemp  ) + &
        fmajor(1,2,1) * k(igpt, jeta(1)  , jpress  , jtemp  ) + &
        fmajor(2,2,1) * k(igpt, jeta(1)+1, jpress  , jtemp  ) ) + &
      scaling(2) * &
      ( fmajor(1,1,2) * k(igpt, jeta(2)  , jpress-1, jtemp+1) + &
        fmajor(2,1,2) * k(igpt, jeta(2)+1, jpress-1, jtemp+1) + &
        fmajor(1,2,2) * k(igpt, jeta(2)  , jpress  , jtemp+1) + &
        fmajor(2,2,2) * k(igpt, jeta(2)+1, jpress  , jtemp+1) )
  end function interpolate3D
  ! ----------------------------------------------------------
  pure function interpolate3D_byflav(scaling, fmajor, k, gptS, gptE, jeta, jtemp, jpress) result(res)
    real(wp), dimension(2),     intent(in) :: scaling
    real(wp), dimension(2,2,2), intent(in) :: fmajor ! interpolation fractions for major species
                                                     ! index(1) : reference eta level (temperature dependent)
                                                     ! index(2) : reference pressure level
                                                     ! index(3) : reference temperature level
    real(wp), dimension(:,:,:,:),intent(in) :: k ! (gpt, eta,temp,press)
    integer,                     intent(in) :: gptS, gptE
    integer, dimension(2),       intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    integer,                     intent(in) :: jtemp ! interpolation index for temperature
    integer,                     intent(in) :: jpress ! interpolation index for pressure
    real(wp), dimension(gptE-gptS+1)        :: res ! the result

    ! Local variable
    integer :: igpt
    ! each code block is for a different reference temperature
    do igpt = 1, gptE-gptS+1
      res(igpt) =  &
        scaling(1) * &
        ( fmajor(1,1,1) * k(gptS+igpt-1, jeta(1)  , jpress-1, jtemp  ) + &
          fmajor(2,1,1) * k(gptS+igpt-1, jeta(1)+1, jpress-1, jtemp  ) + &
          fmajor(1,2,1) * k(gptS+igpt-1, jeta(1)  , jpress  , jtemp  ) + &
          fmajor(2,2,1) * k(gptS+igpt-1, jeta(1)+1, jpress  , jtemp  ) ) + &
        scaling(2) * &
        ( fmajor(1,1,2) * k(gptS+igpt-1, jeta(2)  , jpress-1, jtemp+1) + &
          fmajor(2,1,2) * k(gptS+igpt-1, jeta(2)+1, jpress-1, jtemp+1) + &
          fmajor(1,2,2) * k(gptS+igpt-1, jeta(2)  , jpress  , jtemp+1) + &
          fmajor(2,2,2) * k(gptS+igpt-1, jeta(2)+1, jpress  , jtemp+1) )
    end do
  end function interpolate3D_byflav
  ! ----------------------------------------------------------
  !
  ! Combine absoprtion and Rayleigh optical depths for total tau, ssa, g
  !
  pure subroutine combine_and_reorder_2str(ncol, nlay, ngpt, tau_abs, tau_rayleigh, tau, ssa, g) &
      bind(C, name="combine_and_reorder_2str")
    integer,                             intent(in) :: ncol, nlay, ngpt
    real(wp), dimension(ngpt,nlay,ncol), intent(in   ) :: tau_abs, tau_rayleigh
    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau, ssa, g ! inout because components are allocated
    ! -----------------------
    integer  :: icol, ilay, igpt
    real(wp) :: t
    ! -----------------------
    do icol = 1, ncol
      do ilay = 1, nlay
        do igpt = 1, ngpt
           t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
           tau(icol,ilay,igpt) = t
           g  (icol,ilay,igpt) = 0._wp
           if(t > 2._wp * tiny(t)) then
             ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t
           else
             ssa(icol,ilay,igpt) = 0._wp
           end if
        end do
      end do
    end do
  end subroutine combine_and_reorder_2str
  ! ----------------------------------------------------------
  !
  ! Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
  !   using Rayleigh scattering phase function
  !
  pure subroutine combine_and_reorder_nstr(ncol, nlay, ngpt, nmom, tau_abs, tau_rayleigh, tau, ssa, p) &
      bind(C, name="combine_and_reorder_nstr")
    integer, intent(in) :: ncol, nlay, ngpt, nmom
    real(wp), dimension(ngpt,nlay,ncol), intent(in ) :: tau_abs, tau_rayleigh
    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau, ssa
    real(wp), dimension(ncol,nlay,ngpt,nmom), &
                                         intent(inout) :: p
    ! -----------------------
    integer :: icol, ilay, igpt, imom
    real(wp) :: t
    ! -----------------------
    do icol = 1, ncol
      do ilay = 1, nlay
        do igpt = 1, ngpt
          t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
          tau(icol,ilay,igpt) = t
          if(t > 2._wp * tiny(t)) then
            ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t
          else
            ssa(icol,ilay,igpt) = 0._wp
          end if
          do imom = 1, nmom
            p(imom,icol,ilay,igpt) = 0.0_wp
          end do
          if(nmom >= 2) p(2,icol,ilay,igpt) = 0.1_wp
        end do
      end do
    end do
  end subroutine combine_and_reorder_nstr
end module mo_gas_optics_kernels
