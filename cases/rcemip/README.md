# RCEMIP radiative convective equilibrium case.
This case contains the tropical convection case as documented in the paper of *Wing et al. (2018): Radiative Convective Equilibrium Intercomparison Project (RCEMIP), GMD.*

In order to run this case, it is necessary to copy (or link) the spectral coefficients NetCDF files from the repository of the RTE+RRTMGP radiation (https://github.com/RobertPincus/rte-rrtmgp) code into this directory. Note that you need to use the `v1.0.0` release of RTE+RRTMGP, which can be obtained using:

    git checkout tags/v1.0.0

The gas optics (https://github.com/RobertPincus/rte-rrtmgp/tree/master/rrtmgp/data) data should be linked as `coefficients_lw.nc` and `coefficients_sw.nc` for longwave and shortwave radiation respectivily. Coefficents for cloud optics (https://github.com/RobertPincus/rte-rrtmgp/tree/master/extensions/cloud_optics) should be linked as `cloud_coefficients_lw.nc` and `cloud_coefficients_sw.nc`.
