# EUREC4A model intercomparison (open boundary setup)

https://eurec4a.eu/mip

## Scripts / order.

1. Download required data.
- `download_cosmo.py`: Download 2D/3D COSMO data from swift.dkrz.de.
- `download_era5.py`: Download ERA5 data using (LS)2D, required for RRTMGP background.

2. Pre-process data.
- `create_background.py`: Create time dependent mean background profiles (needed for e.g. RRTMGP), by blending COSMO (below ~20 km) and ERA5 (from ~20 km to TOA). This results in a `eurec4a_mean_profiles.nc` NetCDF file.

3. Domain definitions.
- `grid_definition.py` and `grid_definition_dev.py` define the domains/nesting. The `_dev` version uses a small domain for quick testing.

4. Generate case input.
- `create_init_lbcs_outer.py` creates the initial 3D fields and lateral/top boundary conditions (`thl`, `qt`, `qr`, `u`, `v`, `w`) from COSMO. All variables except `w` are tri-linearly interpolated, after which `w` is calculated to satisfy `∂(ρui)/∂xi ≡ 0`.
- `create_init_lbcs_inner.py` TODO.
- `create_input.py` -- used for both the inner and outer domain -- creates the `eurec4a_input.nc` file needed by MicroHH, containing e.g. the RRTMGP radiation background profiles. Additionally, this script creates the time varying 2D `thl_bot` and `qt_bot` fields, and adds the right settings from `eurec4a.ini.base` to and new `eurec4a.ini` file.

NOTE: `create_init_lbcs_outer.py` and `create_init_lbcs_inner.py` write the initial 3D fields to files named e.g `u_0.0000000`. Between the MicroHH `init` and `run` phase, these files need to be moved to e.g. `u.0000000` to overwrite the 3D files created by the `init` phase.

