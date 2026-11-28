# Power plant plume with chemistry

## Case description
Idealised version of the Jänschwalde (Germany) power plant case from the CoCO2 project. The emissions (CO2, NO, NO2, CO) of the full plant are distributed over the nine cooling towers, and released as Gaussian blobs.

The case has several options, which can be set at the top of `plume_chem_input.py`:

- `sw_plume_rise`: resolved plume rise, by emitting heat and moisture from the cooling towers. Without plume rise, the emissions are released at a fixed height fitted to the CoCO2 emission profiles.
- `sw_chemistry`: non-linear KPP chemistry, with photolysis rates from TUV.
- `sw_land_surface`: interactive land-surface model, with deposition coupled to it. Without it, fixed deposition velocities are used.

The scalars always use open (non-periodic) lateral boundaries.

## Required files
`plume_chem_input.py` requires `microhhpy`:

    pip install microhhpy

and the following files need to be linked or copied to the case directory:

    ln -s ../../python/microhh_tools.py
    ln -s ../../python/lsm_input.py
    ln -s ../../misc/van_genuchten_parameters.nc

`van_genuchten_parameters.nc` is only needed with `sw_land_surface=True`.
