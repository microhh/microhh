#!/bin/bash

cd conservation100
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

cd conservation200
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

cd conservation400
python conservationprof.py
rm -f *.00* conservation.out
../init conservation
../dns conservation
cd ..

