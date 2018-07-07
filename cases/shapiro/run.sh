#!/bin/bash
cd 32x64
$PYTHON_EXEC shapiroprof.py 
rm -f *.0* shapiro.out 
../microhh init shapiro
../microhh run shapiro
cd ..

cd 64x128
$PYTHON_EXEC shapiroprof.py 
rm -f *.0* shapiro.out 
../microhh init shapiro
../microhh run shapiro
cd ..

cd 128x256
$PYTHON_EXEC shapiroprof.py 
rm -f *.0* shapiro.out 
../microhh init shapiro
../microhh run shapiro
cd ..
