#!/bin/bash

basedir=$(pwd)
MICROHH_EXEC=$basedir/"../build_gpu/microhh_single" #Choose here to use the CPU, GPU, MPI, etc version of MicroHH
echo "Test runs at" $(date)
echo "MicroHH Binary:" $MICROHH_EXEC
for dir in */
do
    echo $dir
    cd $dir
    if [ -f run.sh ]; then
      time run.sh
    fi
    cd $basedir
done
