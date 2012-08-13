#!/bin/bash

cd conservation100_2nd
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

cd conservation200_2nd
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

cd conservation400_2nd
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

cd conservation100_4th
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

cd conservation200_4th
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

cd conservation400_4th
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

