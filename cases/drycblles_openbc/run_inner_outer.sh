rm inner/*
rm outer/*

rm *00*
python drycblles_input.py outer
./microhh init drycblles
#./microhh run drycblles
gdb --args ./microhh run drycblles

#python cross_to_nc.py -n 6
#mv *.nc outer/
#mv *00* outer/
#
#rm *00*
#python drycblles_input.py inner
#./microhh init drycblles
#./microhh run drycblles
#python cross_to_nc.py -n 6
#mv *.nc inner/
#mv *00* inner/
