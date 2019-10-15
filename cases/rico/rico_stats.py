import xarray as xr
import matplotlib.pyplot as pl

pl.close('all')

f = xr.open_dataset('cases/rico/rico_default_0000000.nc', decode_times=False) # !! CHANGED FROM rico.default.0000000.nc

# Colors in plot
c_sedi = 'C0'    # Sedimentation
c_auto = 'C1'    # Autoconversion
c_evap = 'C2'    # Evaporation
c_scbr = 'C3'    # Selfcollection and breakup
c_accr = 'C4'    # Accretion

# Time (index_ to plot
time = 15

pl.figure()
pl.subplot(221)
pl.xlabel('dthl/dt (K h-1)')
pl.plot(f['auto_thlt'][time,:]*3600., f['z'], label='Autconversion', color=c_auto)
pl.plot(f['evap_thlt'][time,:]*3600., f['z'], label='Evaporation', color=c_evap)
pl.plot(f['accr_thlt'][time,:]*3600., f['z'], label='Accretion', color=c_accr)
pl.legend()

pl.subplot(222)
pl.xlabel('dqt/dt (g kg-1 h-1)')
pl.plot(f['auto_qtt'][time,:]*3600000., f['z'], label='Autconversion', color=c_auto)
pl.plot(f['evap_qtt'][time,:]*3600000., f['z'], label='Evaporation', color=c_evap)
pl.plot(f['accr_qtt'][time,:]*3600000., f['z'], label='Accretion', color=c_accr)
pl.legend()

pl.subplot(223)
pl.xlabel('dqr/dt (g kg-1 h-1)')
pl.plot(f['auto_qrt'][time,:]*3600000., f['z'], label='Autconversion', color=c_auto)
pl.plot(f['evap_qrt'][time,:]*3600000., f['z'], label='Evaporation', color=c_evap)
pl.plot(f['accr_qrt'][time,:]*3600000., f['z'], label='Accretion', color=c_accr)
pl.plot(f['sed_qrt' ][time,:]*3600000., f['z'], label='Sedimentation', color=c_sedi)
pl.legend()

pl.subplot(224)
pl.xlabel('dnr/dt (m-3 h-1)')
pl.plot(f['auto_nrt'][time,:]*3600000., f['z'], label='Autconversion', color=c_auto)
pl.plot(f['evap_nrt'][time,:]*3600000., f['z'], label='Evaporation', color=c_evap)
pl.plot(f['scbr_nrt'][time,:]*3600000., f['z'], label='Selfcolletion/breakup', color=c_scbr)
pl.plot(f['sed_nrt' ][time,:]*3600000., f['z'], label='Sedimentation', color=c_sedi)
pl.legend(loc=1)
