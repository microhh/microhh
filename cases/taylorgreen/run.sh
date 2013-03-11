#!/bin/bash

# 2nd order scheme
cd taylorgreen16_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen32_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen64_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen128_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen256_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

# 4th order scheme Morinishi
cd taylorgreen16_44
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen32_44
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen64_44
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen128_44
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen256_44
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

# 4th order scheme
cd taylorgreen16_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen32_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen64_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen128_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

cd taylorgreen256_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../microhh init taylorgreen
../microhh run taylorgreen
cd ..

