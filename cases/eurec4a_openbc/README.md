# EUREC4A model intercomparison (open boundary setup)

https://eurec4a.eu/mip

## Scripts / order.

0. Settings.
- `global_settings.py`: some settings shared between scripts are defined here.

1. Download required data.
- `download_cosmo.py`: Download 2D/3D COSMO data from swift.dkrz.de.
- `download_era5.py`: Download ERA5 data using (LS)2D, required for RRTMGP background.

2. Pre-process data.
- `create_background.py`: Create time dependent mean background profiles (needed for e.g. RRTMGP), by blending COSMO (below ~20 km) and ERA5 (from ~20 km to TOA). This is somewhat costly, so do it once for the full EUREC4A period. This results in a `eurec4a_mean_profiles.nc` NetCDF file.

3. Domain definitions.
- `grid_definition.py` define the domains/nesting. The switch between the different domains is set in `global_settings.py`.

4. Generate case input.
- `create_default_input.py` -- used for both the inner and outer domain -- creates the `eurec4a_input.nc` file needed by MicroHH, containing e.g. the RRTMGP radiation background profiles. Additionally, this script creates the time varying 2D `thl_bot` and `qt_bot` fields, and adds the right settings from `eurec4a.ini.base` to and new `eurec4a.ini` file. Finally, it generates the thermodynamic base state.
- TODO.

NOTE: `create_init_lbcs_outer.py`, `create_init_lbcs_inner.py`, and `create_default_input.py` write the initial 3D fields and base state to files named `u_0.0000000` or `thermo_basestate_0.0000000`, et cetera. Between the MicroHH `init` and `run` phase, these files need to be moved to e.g. `u.0000000` to overwrite the files created by the `init` phase.