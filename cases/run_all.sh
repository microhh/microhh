#!/bin/bash

basedir=$(pwd)
export PYTHON_EXEC=python3
export MICROHH_EXEC=$basedir/"../build_gpu/microhh_single" #Choose here to use the CPU, GPU, MPI, etc version of MicroHH
blacklist="moser600/ bomex/"
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
