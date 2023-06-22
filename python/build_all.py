import subprocess
import os

def execute(task):
    subprocess.call(task, shell=True, executable='/bin/bash')

TFs = ['sp', 'dp']
archs = ['cpu', 'cpumpi', 'gpu']
build_type = 'release'

for arch in archs:
    for TF in TFs:
        build_dir = '../build_{}_{}'.format(TF, arch)
        if not os.path.exists(build_dir):
            os.mkdir(build_dir)
        os.chdir(build_dir)

        usecuda = (arch == 'gpu')
        usempi  = (arch == 'cpumpi')
        usesp   = (TF == 'sp')

        cmake_cmd = 'cmake -DUSEMPI={} -DUSECUDA={} -DUSESP={} -DCMAKE_BUILD_TYPE={}  ..'.format(
            usempi, usecuda, usesp, build_type)
        execute(cmake_cmd)
        execute('make -j 4')

        os.chdir('../python')
