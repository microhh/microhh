#!/bin/bash

basedir=$(pwd)
export PYTHON_EXEC=python3
export MICROHH_EXEC=$basedir/"../build_gpu/microhh_single" #Choose here to use the CPU, GPU, MPI, etc version of MicroHH
blacklist="moser600/ eady/ eady_callies/ gabls4s3/ gabls4s3_nbl/ vanheerwarden2016/"
echo "Test runs at" $(date)
echo "MicroHH Binary:" $MICROHH_EXEC
for dir in */
do
    echo $dir
    if [[ $blacklist == *$dir* ]]; then
      echo "..Skipping"
    else
      cd $dir
      if [ -f run.sh ]; then
        time run.sh
      fi
   fi
   cd $basedir
     
done
