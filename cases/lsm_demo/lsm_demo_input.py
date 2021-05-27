import matplotlib.pyplot as pl
import netCDF4 as nc4
import numpy as np
from datetime import datetime
import shutil

# Custom scripts
from tools import Read_afgl, Grid_linear_stretched
from lsm_input import LSM_input
import microhh_tools as mht

# Constants
Rd = 287.04
Rv = 461.5
ep = Rd/Rv

# Switch between single and double precision mode:
TF = np.float64

# Vertical grid used by LES and the background radiation profile:
grid_les = Grid_linear_stretched(ktot=64, dz0=20, alpha=0.03)
grid_rad = Grid_linear_stretched(ktot=96, dz0=50, alpha=0.04)
#grid_les.plot()
#grid_rad.plot()

print('zsize={}'.format(grid_les.zsize))

# Read MicroHH namelist & grid properties:
nl = mht.Read_namelist('lsm_demo.ini')

itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot_soil = nl['land_surface']['ktot']

xsize = nl['grid']['xsize']
ysize = nl['grid']['ysize']

dx = xsize/itot
dy = ysize/jtot

# Read AFGL standard atmosphere (midlatitude summer) as background sounding for RRTMGP:
afgl = Read_afgl('afgl_midlat_s.snd')

# --------------------------------------
# Define initial profiles LES and RRTMGP
# --------------------------------------
# Well mixed layer with linear stratification on top:
zi    = 200.
th_b  = 290
dth   = 3.
dthdz = 5e-3
qt_b  = 8e-3
dqt   = -2e-3
dqtdz = -1e-6

# Mask for ABL and FT:
abl = grid_les.z <= zi
ft  = ~abl

th = np.zeros(grid_les.ktot)
qt = np.zeros(grid_les.ktot)
u  = np.zeros(grid_les.ktot)
ug = np.zeros(grid_les.ktot)

th[abl] = th_b
th[ft ] = th_b + dth + (grid_les.z[ft]-zi)*dthdz

qt[abl] = qt_b
qt[ft ] = qt_b + dqt + (grid_les.z[ft]-zi)*dqtdz

u[:] = 5.
ug[:] = 5.

# Interpolate AFGL ozone onto LES grid:
o3 = np.interp(grid_les.z, afgl.z, afgl.o3)

# Interpolate AFGL soundings onto radiation grid:
afgl_p   = np.interp(grid_rad.z, afgl.z, afgl.p)
afgl_th  = np.interp(grid_rad.z, afgl.z, afgl.thl)
afgl_qt  = np.interp(grid_rad.z, afgl.z, afgl.qt)
afgl_o3  = np.interp(grid_rad.z, afgl.z, afgl.o3)
afgl_h2o = np.interp(grid_rad.z, afgl.z, afgl.h2o)

# Set values background profiles inside LES domain:
abl = grid_rad.z <= zi
ft  = (~abl) & (grid_rad.z < grid_les.zsize)

afgl_th[abl] = th_b
afgl_th[ft ] = th_b + dth + (grid_rad.z[ft]-zi)*dthdz

afgl_qt[abl] = qt_b
afgl_qt[ft ] = qt_b + dqt + (grid_rad.z[ft]-zi)*dqtdz

# Blend AFGL soundings to LES values above LES domain:
th_top_les = th_b + dth + (grid_les.zsize-zi)*dthdz
qt_top_les = qt_b + dqt + (grid_les.zsize-zi)*dqtdz

th_top_afgl = np.interp(grid_les.zsize, afgl.z, afgl.thl)
qt_top_afgl = np.interp(grid_les.zsize, afgl.z, afgl.qt)

dth_top = th_top_les - th_top_afgl
dqt_top = qt_top_les - qt_top_afgl

dz_blend = 2000
mask = (grid_rad.z >= grid_les.zsize) & (grid_rad.z <= grid_les.zsize+dz_blend)
fac  = (dz_blend - (grid_rad.z[mask] - grid_les.zsize)) / dz_blend

afgl_th[mask] += fac * dth_top
afgl_qt[mask] += fac * dqt_top

# Calculate background values on half levels:
def calc_half_levels(var_in, z, zsize):
    var_out = np.zeros(var_in.size+1)
    var_out[1:-1] = 0.5 * (var_in[1:] + var_in[:-1])
    var_out[0 ] = var_in[0 ] - (var_in[1 ]-var_in[0 ]) / (z[1 ]-z[0 ]) * z[0]
    var_out[-1] = var_in[-1] + (var_in[-1]-var_in[-2]) / (z[-1]-z[-2]) * (zsize - z[-1])
    return var_out

afgl_ph  = np.interp(grid_rad.zh, afgl.z, afgl.p)
afgl_thh = calc_half_levels(afgl_th, grid_rad.z, grid_rad.zsize)

# Convert to units used by radiation scheme:
Mw = 18.016; Md = 28.966; Rd=287.04; cp=1005; p0=1e5
afgl_T   = afgl_th *  (afgl_p  / p0)**(Rd/cp)
afgl_Th  = afgl_thh * (afgl_ph / p0)**(Rd/cp)
afgl_h2o = afgl_qt / (ep - ep*afgl_qt)

# Soil/surface:
# Midpoint of soil layers:
z_soil    = np.array([-1.945, -0.64, -0.175, -0.035])

# Lookup table with soil properties:
nc_vg = nc4.Dataset('van_genuchten_parameters.nc')  # available in microhh_root/misc

si = 2  # 2=ECMWF medium fine
theta_wp = nc_vg.variables['theta_wp'][si]
theta_fc = nc_vg.variables['theta_fc'][si]

index_soil = np.zeros(ktot_soil)
theta_soil = np.zeros(ktot_soil)
t_soil     = np.zeros(ktot_soil)

index_soil[:] = si
theta_soil[:] = theta_wp + 0.75 * (theta_fc - theta_wp)
t_soil[:] = 291
root_frac = np.array([0.04, 0.23, 0.38, 0.35])

# -------------------------
# Write LES input to NetCDF
# -------------------------
def add_nc_var(nc_group, name, dims, values, dtype=TF):
    var = nc_group.createVariable(name, dtype, dims)
    var[:] = values

nc = nc4.Dataset('lsm_demo_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

# Create groups:
nc_init = nc.createGroup('init')
nc_rad  = nc.createGroup('radiation')
nc_soil = nc.createGroup('soil')

# Create dimensions:
nc.createDimension('z', grid_les.ktot)
nc_rad.createDimension('lay', grid_rad.ktot)
nc_rad.createDimension('lev', grid_rad.ktot+1)
nc_soil.createDimension('z', ktot_soil)

# Write vertical grid:
add_nc_var(nc, 'z', ('z'), grid_les.z)

# Initial profiles:
add_nc_var(nc_init, 'thl', ('z'), th)
add_nc_var(nc_init, 'qt',  ('z'), qt)
add_nc_var(nc_init, 'u',   ('z'), u)
add_nc_var(nc_init, 'ug',  ('z'), ug)

add_nc_var(nc_init, 'h2o', ('z'), qt/(ep-ep*qt))
add_nc_var(nc_init, 'o3',  ('z'), o3*1e-6)
add_nc_var(nc_init, 'co2', (), 348e-6)
add_nc_var(nc_init, 'ch4', (), 1650e-9)
add_nc_var(nc_init, 'n2o', (), 306e-9)
add_nc_var(nc_init, 'n2',  (), 0.7808)
add_nc_var(nc_init, 'o2',  (), 0.2095)

# Radiation profiles:
nc_rad = nc.createGroup('radiation')
add_nc_var(nc_rad, 'z_lay', ('lay'), grid_rad.z)
add_nc_var(nc_rad, 'z_lev', ('lev'), grid_rad.zh)

add_nc_var(nc_rad, 'p_lay', ('lay'), afgl_p)
add_nc_var(nc_rad, 'p_lev', ('lev'), afgl_ph)

add_nc_var(nc_rad, 't_lay', ('lay'), afgl_T)
add_nc_var(nc_rad, 't_lev', ('lev'), afgl_Th)

add_nc_var(nc_rad, 'h2o', ('lay'), afgl_h2o)
add_nc_var(nc_rad, 'o3',  ('lay'), afgl_o3*1e-6)
add_nc_var(nc_rad, 'co2', (), 348e-6)
add_nc_var(nc_rad, 'ch4', (), 1650e-9)
add_nc_var(nc_rad, 'n2o', (), 306e-9)
add_nc_var(nc_rad, 'n2',  (), 0.7808)
add_nc_var(nc_rad, 'o2',  (), 0.2095)

# Soil profiles:
add_nc_var(nc_soil, 'z', ('z'), z_soil)

if nl['land_surface']['swhomogeneous']:
    add_nc_var(nc_soil, 't_soil', ('z'), t_soil)
    add_nc_var(nc_soil, 'theta_soil', ('z'), theta_soil)
    add_nc_var(nc_soil, 'root_frac', ('z'), root_frac)
    add_nc_var(nc_soil, 'index_soil', ('z'), index_soil, np.int32)

nc.close()

# -------------------------------------------
# Create input for heterogeneous land surface
# -------------------------------------------
if not nl['land_surface']['swhomogeneous']:

    # Help class to simplify writing the input files for heterogeneous land surfaces:
    lsm_input = LSM_input(
            itot, jtot, nl['land_surface']['ktot'], TF=TF, debug=True)

    # Initialise all values:
    lsm_input.c_veg[:,:] = 0.95
    lsm_input.z0m[:,:] = 0.075
    lsm_input.z0h[:,:] = 0.003
    lsm_input.gD[:,:] = 0.
    lsm_input.lai[:,:] = 2.6
    lsm_input.rs_veg_min[:,:] = 100.
    lsm_input.rs_soil_min[:,:] = 50.
    lsm_input.lambda_stable[:,:] = 5.
    lsm_input.lambda_unstable[:,:] = 5.
    lsm_input.cs_veg[:,:] = 0.
    lsm_input.water_mask[:,:] = 0
    lsm_input.t_soil[:,:,:] = t_soil[:,np.newaxis,np.newaxis]
    lsm_input.index_soil[:,:,:] = index_soil[:,np.newaxis,np.newaxis]
    lsm_input.root_frac[:,:,:] = root_frac[:,np.newaxis,np.newaxis]

    # Create irrigation circles as example..
    n = 4
    D = xsize / n

    # LES spatial grid:
    x = np.arange(dx/2, xsize, dx)
    y = np.arange(dy/2, ysize, dy)
    x2,y2 = np.meshgrid(x,y)

    theta_moist = theta_wp + 0.8 * (theta_fc - theta_wp)
    theta_dry   = theta_wp + 0.3 * (theta_fc - theta_wp)

    lsm_input.theta_soil[:,:,:] = theta_dry

    for i in range(n):
        for j in range(n):
            x0 = (i+0.5)*D
            y0 = (j+0.5)*D
            d = np.sqrt((x2-x0)**2 + (y2-y0)**2)
            m = d < (D/2)-25.
            lsm_input.theta_soil[:,m] = theta_moist

    # Check if all values are initialised:
    lsm_input.check()

    # Write binary input files for model:
    lsm_input.save_binaries(allow_overwrite=True)

    # Write lsm input as NetCDF file (visualisation, ...):
    lsm_input.save_netcdf('lsm_input.nc', allow_overwrite=True)
