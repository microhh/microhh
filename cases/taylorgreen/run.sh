#!/bin/bash

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

