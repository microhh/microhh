#!/bin/bash

cd ekman32
python ekmanprof.py
rm -f *.0* ekman.out
../microhh init ekman
../microhh run ekman
cd ..

cd ekman64
python ekmanprof.py
rm -f *.0* ekman.out
../microhh init ekman
../microhh run ekman
cd ..

cd ekman128
python ekmanprof.py
rm -f *.0* ekman.out
../microhh init ekman
../microhh run ekman
cd ..

cd ekman256
python ekmanprof.py
rm -f *.0* ekman.out
../microhh init ekman
../microhh run ekman
cd ..

