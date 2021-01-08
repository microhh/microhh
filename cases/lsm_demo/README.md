# Interactive land surface example

Idealised example case with homogeneous or heterogeneous interactive land surface.

1. In order to run this case, it is necessary to copy (or link) the spectral coefficients NetCDF files from the repository of the RTE+RRTMGP radiation (https://github.com/RobertPincus/rte-rrtmgp, the v1.0 release) code into this directory. The gas optics (https://github.com/RobertPincus/rte-rrtmgp/tree/master/rrtmgp/data) data should be linked as `coefficients_lw.nc` and `coefficients_sw.nc` for longwave and shortwave radiation respectivily. Coefficents for cloud optics (https://github.com/RobertPincus/rte-rrtmgp/tree/master/extensions/cloud_optics) should be linked as `cloud_coefficients_lw.nc` and `cloud_coefficients_sw.nc`.
2. Copy (or link) the van Genuchten lookup table (`<microhh_root>/misc/van_genuchten_parameters.nc`) to this directory. This lookup table is used both by the initialisation script, and the model itself.
3. Run `lsm_demo_input.py` to generate the input for MicroHH.

The case is by default configured for a homogeneous land surface, controlled through the namelist options in the `land_surface` group. With `sw_homogeneous=0`, the `lsm_demo_input.py` script generates the binary files needed to initialise the heterogeneous land surface.
