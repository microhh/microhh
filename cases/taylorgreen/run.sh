#!/bin/bash

cd taylorgreen32
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen64
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen128
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen256
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

