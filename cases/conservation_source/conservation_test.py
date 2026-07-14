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


# Point source definition, shared by the `gaussian` and `3d` source methods
# so both test the exact same physical source.
source_x0 = 1200
source_y0 = 800
source_z0 = 200
sigma_x = 50
sigma_y = 50
sigma_z = 50
strength = 1
sw_vmr = False


def test_conservation(sw_thermo, sw_basestate, sw_source, executable, precision):
    """
    Test scalar conservation from point source.
    """
    print(f'Testing {executable} for swthermo={sw_thermo}, swbasestate={sw_basestate}, swsource={sw_source}')
    float_type = np.float32 if precision == 'sp' else np.float64

    """
    1. Update .ini file
    """
    ini = mht.Read_namelist('conservation.ini.base')

    ini['thermo']['swthermo'] = sw_thermo
    ini['thermo']['swbasestate'] = sw_basestate

    thermo_var = 'th' if sw_thermo == 'dry' else 'thl'
    ini['fields'][f'rndamp[{thermo_var}]'] = 0.1

    """
    2. Define grid and input profiles.
    """
    itot = ini['grid']['itot']
    jtot = ini['grid']['jtot']
    ktot = ini['grid']['ktot']

    xsize = ini['grid']['xsize']
    ysize = ini['grid']['ysize']
    zsize = ini['grid']['zsize']

    pbot = ini['thermo']['pbot']

    dx = xsize / itot
    dy = ysize / jtot
    dz = zsize / ktot

    x = np.arange(dx/2, xsize, dx)
    y = np.arange(dy/2, ysize, dy)
    z = np.arange(dz/2, zsize, dz)

    s = np.zeros(ktot)
    u = np.ones(ktot)*5
    th = 290 + z * 0.006
    qt = np.zeros(ktot)

    """
    3. Configure source, and generate 3D emission binaries if needed.
    """
    if sw_source == 'gaussian':
        ini['source']['swsource'] = 'gaussian'
        ini['source']['sourcelist'] = 's1'
        ini['source']['source_x0'] = source_x0
        ini['source']['source_y0'] = source_y0
        ini['source']['source_z0'] = source_z0
        ini['source']['sigma_x'] = sigma_x
        ini['source']['sigma_y'] = sigma_y
        ini['source']['sigma_z'] = sigma_z
        ini['source']['strength'] = strength
        ini['source']['swvmr'] = sw_vmr
        emiss = None

    elif sw_source == '3d':
        # Only `import microhhpy` for the 3D source method, so the other
        # tests can run without it being installed.
        from microhhpy.thermo import calc_dry_basestate, calc_moist_basestate
        from microhhpy.chem import Emission_input

        # Base state density, matching MicroHH exactly. `Emission_input`
        # normalizes the emission with this density, so it has to be correct.
        if sw_basestate == 'boussinesq':
            rho_ref = np.ones(ktot, float_type)
        elif sw_thermo == 'dry':
            rho_ref = calc_dry_basestate(th, pbot, z, zsize, float_type)['rho']
        else:
            rho_ref = calc_moist_basestate(th, qt, pbot, z, zsize, float_type)['rho']

        # Create emission input instance, and add the Gaussian source.
        emiss = Emission_input(['s1'], [0], x, y, z, np.full(ktot, dz), rho_ref, float_type)

        emiss.add_gaussian(
                's1', strength, 0, source_x0, source_y0, source_z0,
                sigma_x, sigma_y, sigma_z, sw_vmr)

        # Clip to required vertical extent.
        emiss.clip()

        ini['source']['swsource'] = '3d'
        ini['source']['sourcelist'] = 's1'
        ini['source']['ktot'] = emiss.kmax
        ini['source']['strength'] = strength

    ini.save('conservation.ini', allow_overwrite=True)

    """
    4. Create NetCDF case input.
    """
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
    5. Run case
    """
    clean_case()   # Just to be sure case can start.

    # Write emission binaries after `clean_case()`, as they match `*00*`.
    if emiss is not None:
        emiss.to_binary(path='.')

    status = 0
    status += execute(f'{executable} init conservation')
    status += execute(f'{executable} run conservation')

    if status > 0:
        print(f'Running case with executable {executable} failed!')
    else:

        """
        6. Check mass conservation for periodic scalar.
        """
        endtime = ini['time']['endtime']
        savetime = ini['time']['savetime']

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


def run_conservation_test(modes, precs, sources, thermos, bases):

    status = 0

    for sw_source in sources:
        print('-----------------')
        print(f'Source mode = {sw_source}')
        print('-----------------')

        for mode in modes:
            for prec in precs:
                executable = f'../../build_{prec}_{mode}/microhh'
                for sw_thermo in thermos:
                    for sw_basestate in bases:
                        status += test_conservation(
                                sw_thermo, sw_basestate, sw_source, executable, prec)

    return status


if __name__ == '__main__':
    """
    Run full conservation test.
    """

    modes = ['cpu', 'cpumpi', 'gpu']
    precs = ['sp', 'dp']
    sources = ['gaussian', '3d']
    thermos = ['dry', 'moist']
    bases = ['boussinesq', 'anelastic']

    status = run_conservation_test(modes, precs, sources, thermos, bases)

    if status > 0:
        raise Exception('One or more conservation tests failed.')
