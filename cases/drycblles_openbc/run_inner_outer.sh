rm inner/*
rm outer/*

rm *00*
python drycblles_input.py outer
mpiexec -n 8 ./microhh init drycblles
mpiexec -n 8 ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc outer/
mv *00* outer/

rm *00*
python drycblles_input.py inner
mpiexec -n 8 ./microhh init drycblles
mpiexec -n 8 ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc inner/
mv *00* outer/
