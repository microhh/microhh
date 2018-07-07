#!/bin/bash

cd conservation100_3rd
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation200_3rd
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation400_3rd
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation800_3rd
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation100_4th
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation200_4th
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation400_4th
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation800_4th
$PYTHON_EXEC conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..
