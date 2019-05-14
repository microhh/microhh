import numpy as np
import netCDF4 as nc

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open('rcemip.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

thl = np.zeros(np.size(z))
qt  = np.zeros(np.size(z))

nc_file = nc.Dataset("rcemip_input.nc", mode="w", datamodel="NETCDF4", clobber=False)

# The height dimension.
nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))
nc_z[:] = z[:]

# Initial profiles.
nc_group_init = nc_file.createGroup("init");
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt  = nc_group_init.createVariable("qt" , float_type, ("z"))

nc_thl[:] = thl[:]
nc_qt [:] = qt [:]

# Radiation profiles.
p0 = 101480.
p = np.linspace(p0, 0., 1001)

co2 =  348.e-6
ch4 = 1650.e-9
n2o =  306.e-9

g1 = 3.6478
g2 = 0.83209
g3 = 11.3515
p_hpa = p/100.
o3 = g1 * p_hpa**g2 * np.exp(-p_hpa/g3) * 1e-6

nc_file.createDimension("p", p.size)
nc_p = nc_file.createVariable("p", float_type, ("p"))
nc_p[:] = p[:]

nc_group_rad = nc_file.createGroup("radiation")

nc_CO2 = nc_group_rad.createVariable("co2", float_type, ("p"))
nc_CH4 = nc_group_rad.createVariable("ch4", float_type, ("p"))
nc_N2O = nc_group_rad.createVariable("n2o", float_type, ("p"))
nc_O3  = nc_group_rad.createVariable("o3" , float_type, ("p"))

nc_CO2[:] = co2
nc_CH4[:] = ch4
nc_N2O[:] = n2o
nc_O3 [:] = o3[:]

nc_file.close()
