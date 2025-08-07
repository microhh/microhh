import numpy as np
import subprocess
import shutil
import glob
import os

from microhhpy.io import read_ini, save_ini


def is_multiple_of(number, x):
    if x == 0:
        return True
    return number % x == 0


def setup_time():
    """
    Setup start and end of simulation period.
    """

    files = glob.glob('time.0*')
    files.sort()
    restart_time = int(files[-1].split('.')[-1])

    if not is_multiple_of(restart_time, time_chunk):
        raise ValueError('Oink, thats not good.')

    chunck = 0 if restart_time == 0 else restart_time // time_chunk

    start_time_sim = chunck * time_chunk
    end_time_sim = start_time_sim + time_chunk

    return chunck, start_time_sim, end_time_sim


def update_ini(start_time, end_time):
    """
    Update .ini file with correct start/end time.
    """

    ini = read_ini('eurec4a.ini.base')

    ini['time']['starttime'] = start_time
    ini['time']['endtime'] = end_time

    save_ini(ini, 'eurec4a.ini')


def create_runscript(chunck, ntasks, wc_limit):

    script_name = f'run_{chunck}.slurm'

    with open(script_name, 'w') as f:

        f.write('#!/bin/bash\n')
        f.write('set -e\n\n')

        f.write(f'#SBATCH --job-name=eu4a_o\n')
        f.write(f'#SBATCH --output=mhh-%j.out\n')
        f.write(f'#SBATCH --error=mhh-%j.err\n')
        f.write(f'#SBATCH --partition=genoa\n')
        f.write(f'#SBATCH -n {ntasks}\n')
        f.write(f'#SBATCH --cpus-per-task=1\n')
        f.write(f'#SBATCH --ntasks-per-core=1\n')
        f.write(f'#SBATCH -t {wc_limit}\n')
        f.write(f'#SBATCH --mail-user=bart.vanstratum@wur.nl\n') 
        f.write(f'#SBATCH --mail-type=FAIL\n\n') 

        
        f.write('source ~/setup_env.sh\n\n')
        
        f.write('export OMPI_MCA_fcoll="two_phase"\n')
        f.write('export OMPI_MCA_io_ompio_bytes_per_agg="512MB"\n\n')
        
        if chunck == 0:
            f.write('srun ./microhh init eurec4a\n\n')

            f.write("find . -maxdepth 1 -type f -name '*_ext*' | while read -r file; do\n")
            f.write('    newname="${file/_ext/}"\n')
            f.write('    echo "Renaming: $file -> $newname"\n')
            f.write('    mv "$file" "$newname"\n')
            f.write('done\n\n')
            
        f.write('srun ./microhh run eurec4a\n\n')

        f.write('python autosubmit_case.py\n')

    return script_name


if __name__ == '__main__':

    """
    Global settings.
    """
    start_time = 0
    end_time = 12 * 24 * 3600
    time_chunk = 43200
    ntasks = 6144
    wc_limit = '72:00:00'

    # Make backup of .ini once.
    if not os.path.exists('eurec4a.ini.base'):
        shutil.copyfile('eurec4a.ini', 'eurec4a.ini.base')

    # Find start- and entime of simulation.
    start_time_chunck, end_time_chunck = setup_time()
    
    # Update ini.
    update_ini(start_time_chunck, end_time_chunck)
    
    # Create SLURM run script.
    slurm_script = create_runscript(1, ntasks, wc_limit)

    # Submit new run!
    #subprocess.run('sbatch', slurm_script)
