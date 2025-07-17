rm inner/*
rm outer/*

rm *00*
python drycblles_input.py outer
mpiexec -n 4 ./microhh init drycblles
mpiexec -n 4 ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc outer/
mv *00* outer/

python drycblles_input.py inner
cp outer/lbc_* .
mpiexec -n 4 ./microhh init drycblles
mpiexec -n 4 ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc inner/
mv *00* inner/
