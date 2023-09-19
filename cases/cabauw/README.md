# Shallow cumulus over Cabauw

## Case description
This case features a day with shallow cumulus over Cabauw (the Netherlands) on 15-08-2016, running from 06:00 to 18:00 UTC. The initial conditions and large-scale forcings are from ERA5, obtained using the (LS)<sup>2</sup>D Python package (van Stratum et al. 2023: [10.22541/essoar.168167362.25141062/v1](10.22541/essoar.168167362.25141062/v1)).

This directory contains two cases:
- The `cabauw` case, which is a simplified version of one of the cases from Tijhuis et al. (2022): [https://doi.org/10.1029/2022MS003262](https://doi.org/10.1029/2022MS003262).
- The `cabauw_rt` case, from Veerman et al. (2022): [https://doi.org/10.1029/2022GL100808](https://doi.org/10.1029/2022GL100808).

## Running the cases
Running either one of the cases requires the `microhh_tools.py` help script from the `microhh_root/python` directory.

Running `cabauw_input.py` will generate the `cabauw_input.nc` and `cabauw.ini` input files for MicroHH.

The `cabauw_rt_input.py` script will generate all the cases and sensitivity experiments from Veerman et al. (2022) in separate directories.
