#!/bin/bash
cd 32x64
python shapiroprof.py 
rm -f *.0* shapiro.out 
../microhh init shapiro
../microhh run shapiro
cd ..

cd 64x128
python shapiroprof.py 
rm -f *.0* shapiro.out 
../microhh init shapiro
../microhh run shapiro
cd ..
