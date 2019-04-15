#!/bin/bash
cd 32x64
$PYTHON_EXEC shapiroprof.py 
rm -f *.0* shapiro.out 
$MICROHH_EXEC init shapiro
$MICROHH_EXEC run shapiro
cd ..

cd 64x128
$PYTHON_EXEC shapiroprof.py 
rm -f *.0* shapiro.out 
$MICROHH_EXEC init shapiro
$MICROHH_EXEC run shapiro
cd ..

cd 128x256
$PYTHON_EXEC shapiroprof.py 
rm -f *.0* shapiro.out 
$MICROHH_EXEC init shapiro
$MICROHH_EXEC run shapiro
cd ..
