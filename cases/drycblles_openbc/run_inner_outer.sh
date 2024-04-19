np=4

rm *00*
python drycblles_input.py outer
mpiexec -n $np ./microhh init drycblles
mpiexec -n $np ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc outer/

rm *00*
python drycblles_input.py inner
mpiexec -n $np ./microhh init drycblles
mpiexec -n $np ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc inner/
