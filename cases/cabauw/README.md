# Shallow cumulus over Cabauw

## Case description
This case features a day with shallow cumulus over Cabauw (the Netherlands) on 15-08-2016, running from 06:00 to 18:00 UTC. The initial conditions and large-scale forcings are from ERA5, obtained using the (LS)<sup>2</sup>D Python package: *TODO: add reference..*

This directory contains two cases:
- The `cabauw` case, which is a simplified version of one of the cases from Tijhuis et al. (2022): [dx.doi.org/10.1002/essoar.10511758.1](dx.doi.org/10.1002/essoar.10511758.1)
- The `cabauw_rt` case, from Veerman et al. (2022): *TODO: add reference once pre-print is online..*

## Running the cases
Running either one of the cases requires the `microhh_tools.py` help script from the `microhh_root/python` directory.

Running `cabauw_input.py` will generate the `cabauw_input.nc` and `cabauw.ini` input files for MicroHH. This case currently requires the [develop](https://github.com/microhh/microhh/tree/develop) MicroHH branch.

The `cabauw_rt_input.py` script will generate all the cases and sensitivity experiments from Veerman et al. (2022) in separate directories. This case currently requires the [develop_rt](https://github.com/microhh/microhh/tree/develop_rt) MicroHH branch for the ray tracer.
