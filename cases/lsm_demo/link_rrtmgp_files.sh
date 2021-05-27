path=/home/bart/meteo/models/rte-rrtmgp-cpp/rte-rrtmgp

cp $path/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc coefficients_lw.nc
cp $path/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc coefficients_sw.nc
cp $path/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc  cloud_coefficients_lw.nc
cp $path/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc  cloud_coefficients_sw.nc
