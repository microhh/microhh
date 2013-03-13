#!/bin/bash

cd conservation100_3rd
python conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation200_3rd
python conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation400_3rd
python conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation100_4th
python conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation200_4th
python conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

cd conservation400_4th
python conservationprof.py
rm -f *.00* conservation.out
../microhh init conservation
../microhh run conservation
cd ..

