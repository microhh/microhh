Example case with time varying (surface and large scale) forcings

`[force]`:
- Large scale advection of moisture (`swls` & `swtimedep_ls`)
- Nudging of a scalar "s" (`swnudge` & `swtimedep_nudge`)
- Subsidence (`swwls` & `swtimedep_wls`)
- Time varying geostrophic wind (`swtimedep_geo`)

`[thermo]`:
- Surface pressure (`swtimedep_pbot`)

`[boundary]`:
- Time varying surface fluxes of thl (`swtimedep`)
