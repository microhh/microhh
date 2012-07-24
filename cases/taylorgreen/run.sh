#!/bin/bash

# 2nd order scheme
cd taylorgreen32_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen64_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen128_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen256_2nd
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

# 42 order scheme
cd taylorgreen32_42
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen64_42
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen128_42
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen256_42
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

# 4th order scheme
cd taylorgreen32_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen64_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen128_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen256_4th
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

