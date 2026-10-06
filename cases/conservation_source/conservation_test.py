import numpy as np
import netCDF4 as nc
import xarray as xr
import subprocess
import glob
import os

import microhh_tools as mht

from microhh_tools import execute
#def execute(task):
#    return subprocess.call(task, shell=True, executable='/bin/bash')

def xr_read_all(f, groups=['default'], decode_times=False):
    # Read all NetCDF groups into a single Dataset.
    dss = [xr.open_dataset(f, decode_times=decode_times)]
    for group in groups:
        dss.append(xr.open_dataset(f, group=group, decode_times=decode_times))
    return xr.merge(dss)


def clean_case():
    """
    Clean working directory.
    """
    files = glob.glob('*00*')
    for f in files:
        os.remove(f)


def test_conservation(sw_thermo, sw_basestate, executable, precision):
    """
    Test scalar conservation from point source.
    """
    print(f'Testing {executable} for swthermo={sw_thermo} with swbasestate={sw_basestate}')
    float_type = np.float32 if precision == 'sp' else np.float64

    """
    1. Update .ini file
    """
    ini = mht.Read_namelist('conservation.ini.base')

    ini['thermo']['swthermo'] = sw_thermo
    ini['thermo']['swbasestate'] = sw_basestate

    thermo_var = 'th' if sw_thermo == 'dry' else 'thl'
    ini['fields'][f'rndamp[{thermo_var}]'] = 0.1

    ini.save('conservation.ini', allow_overwrite=True)


    """
    2. Create NetCDF case input.
    """
    itot = ini['grid']['itot']
    jtot = ini['grid']['jtot']
    ktot = ini['grid']['ktot']

    xsize = ini['grid']['xsize']
    ysize = ini['grid']['ysize']
    zsize = ini['grid']['zsize']

    dx = xsize / itot
    dy = ysize / jtot
    dz = zsize / ktot

    z = np.arange(dz/2, zsize, dz)

    s = np.zeros(ktot)
    u = np.ones(ktot)*5
    th = 290 + z * 0.006
    qt = np.zeros(ktot)

    def add_var(name, dims, values, nc_group):
        nc_var = nc_group.createVariable(name, float_type, dims)
        nc_var[:] = values
    
    nc_file = nc.Dataset('conservation_input.nc', mode='w', datamodel='NETCDF4')
    nc_file.createDimension('z', ktot)
    add_var('z',  ('z'), z,  nc_file)
    
    nc_init = nc_file.createGroup('init');
    add_var('u',  ('z'), u,  nc_init)

    if sw_thermo == 'dry':
        add_var('th', ('z'), th, nc_init)
    elif sw_thermo == 'moist':
        add_var('thl', ('z'), th, nc_init)
        add_var('qt', ('z'), qt, nc_init)

    nc_file.close()


    """
    3. Run case
    """
    clean_case()   # Just to be sure case can start.

    status = 0
    status += execute(f'{executable} init conservation')
    status += execute(f'{executable} run conservation')

    if status > 0:
        print(f'Running case with executable {executable} failed!')
    else:

        """
        4. Check mass conservation for periodic scalar.
        """
        endtime = ini['time']['endtime']
        savetime = ini['time']['savetime']
        strength = ini['source']['strength']

        ds = xr_read_all('conservation.default.0000000.nc')

        rhoref = ds.rhoref.values
        zh = ds.zh.values
        dz = zh[1:] - zh[:-1]

        ret = 0
        for time in range(0, endtime+1, savetime):
            expected_mass = strength * time

            fld = np.fromfile(f's1.{time:07d}', float_type)
            fld = fld.reshape((ktot, jtot, itot))
            mass = np.sum(rhoref[:,None,None] * fld * dx * dy * dz[:,None,None])

            if not np.isclose(expected_mass, mass, rtol=1e-5):
                print(f'Mass not conserved! Expected={expected_mass} kg, integral field={mass} kg.')
                status += 1

        clean_case()   # Don't leave messy case behind if last test.

    return status


def run_conservation_test(modes, precs):

    sw_thermos = ['dry', 'moist']
    sw_basestates = ['boussinesq', 'anelastic']

    status = 0

    for mode in modes:
        for prec in precs:
            executable = f'../../build_{prec}_{mode}/microhh'
            
            for sw_thermo in sw_thermos:
                for sw_basestate in sw_basestates:
                    status += test_conservation(sw_thermo, sw_basestate, executable, prec)

    return status


if __name__ == '__main__':
    """
    Run full conservation test including GPU.
    """

    status = run_conservation_test(modes=['cpu', 'gpu'], precs=['sp', 'dp'])

    if status > 0:
        raise Exception('One or more conservation tests failed.')