#!/bin/bash

# 2nd order scheme
cd taylorgreen16_2nd
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen32_2nd
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen64_2nd
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen128_2nd
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen256_2nd
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

# 4th order scheme Morinishi
cd taylorgreen16_4m
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen32_4m
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen64_4m
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen128_4m
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen256_4m
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

# 4th order scheme
cd taylorgreen16_4th
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen32_4th
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen64_4th
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen128_4th
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen256_4th
$PYTHON_EXEC taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

