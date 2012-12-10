#!/bin/bash

# 4th order scheme
cd taylorgreen100
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen200
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

cd taylorgreen400
python taylorgreenprof.py
rm -f *.00* taylorgreen.out
../init taylorgreen
../dns taylorgreen
cd ..

