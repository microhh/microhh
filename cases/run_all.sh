#!/bin/bash

MICROHH_EXEC="../build_gpu/microhh_single" #Choose here to use the CPU, GPU, MPI, etc version of MicroHH
basedir=$(pwd)

for dir in */
do
    echo $dir
    cd $dir
    if [ -f run.sh ]; then
      run.sh
    fi
    cd $basedir
done
