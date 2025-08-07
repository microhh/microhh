import subprocess
import argparse
from pathlib import Path
import glob
import shutil

def run_cross_to_nc(
        t0,
        t1,
        tstep,
        x_values,
        variables,
        dest_dir):
    print(f'Running cross_to_nc.py with tstep={tstep}, x={x_values}, var={variables}...')
    
    cmd = [
        'python', 'cross_to_nc.py', 
        '-n', '16', 
        '-m', 'xy', 
        '-p', 'single',
        '-t0', str(t0), 
        '-t1', str(t1), 
        '-tstep', str(tstep)
    ]
    
    cmd += ['-x'] + [str(x) for x in x_values] + ['-v'] + variables
    
    subprocess.run(cmd, check=True)
    
    dest_path = Path(dest_dir)
    dest_path.mkdir(parents=True, exist_ok=True)
    
    for var in variables:
        src_file = Path(f'{var}.xy.nc')
        dst_file = dest_path / f'{var}.xy.nc'
        
        if src_file.exists():
            src_file.rename(dst_file)
            #print(f'Moved {src_file} to {dst_file}')
        else:
            print(f'Warning: {src_file} not found')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Post-process cross-sections.')
    parser.add_argument('-t0', type=int, required=True, help='Start time')
    parser.add_argument('-t1', type=int, required=True, help='End time')
    args = parser.parse_args()

    dest_dir = f'/gpfs/work2/0/nwo21036/bart/eurec4a_full/output/dom0/{args.t0:07d}'


    """
    Copy statistics.
    """
    stats = [f'eurec4a.default.{args.t0:07d}.nc'] + glob.glob(f'eurec4a.column.*.{args.t0:07d}.nc')
    dest_path = Path(dest_dir) / f'statistics'
    dest_path.mkdir(parents=True, exist_ok=True)

    for f in stats:
        src_file = Path(f)
        dst_file = dest_path / f

        if src_file.exists():
            shutil.copy2(src_file, dst_file)
            #src_file.rename(dst_file)
            print(f'Copied {src_file} to {dst_file}')
        else:
            print(f'Warning: {src_file} not found')



    """
    Post-process and copy cross-sections.
    """
    variable_sets = [
        # TABLE 3:
        {
            'tstep': 3600,
            'x_values': [1, 128],   # NOTE: 0 in dom1
            'variables': [
                'sw_flux_dn',
                'sw_flux_up',
                'lw_flux_dn',
                'lw_flux_up',
                'sw_flux_dn_clear',
                'sw_flux_up_clear',
                'lw_flux_dn_clear',
                'lw_flux_up_clear'],
            'dest_dir': f'{dest_dir}/cross_sfc_tod_3600'
        },
        {
            'tstep': 3600,
            'x_values': [0],
            'variables': [
                'p_bot',
                'thl_bot',
                'thl_fluxbot',
                'qt_fluxbot',
                'u_fluxbot',
                'v_fluxbot'],
            'dest_dir': f'{dest_dir}/cross_sfc_tod_3600'
        },
        # TABLE 4:
        {
            'tstep': 300,
            'x_values': [0],
            'variables': [
                'qt_path',
                'ql_path',
                'qi_path',
                'qr_path',
                'cape',
                'u',
                'v',
                'thl',
                'qt',
                'zi',
                'rr_bot',
                'p_bot'],   # Bonus, for converting thl->T.
            'dest_dir': f'{dest_dir}/cross_10m_path_300'
        },
        # TABLE 5:
        {
            'tstep': 3600,
            'x_values': [18, 40, 56, 91, 118],
            'variables': [
                'u',
                'v',
                'thl',
                'qt'],
            'dest_dir': f'{dest_dir}/cross_z_3600'
        },
        {
            'tstep': 3600,
            'x_values': [18, 41, 56, 92, 118],
            'variables': ['w'],
            'dest_dir': f'{dest_dir}/cross_z_3600'
        },
        {
            'tstep': 3600,
            'x_values': [18, 40, 56],
            'variables': [
                'ql',
                'qr',
                'rh'],
            'dest_dir': f'{dest_dir}/cross_z_3600'
        }
    ]
    
    for var_set in variable_sets:
        run_cross_to_nc(args.t0, args.t1, **var_set)
