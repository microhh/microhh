rm inner/*
rm outer/*

# Outer domain.
rm *00*
python drycblles_input.py outer
mpiexec -n 8 ./microhh init drycblles
mpiexec -n 8 ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc outer/
mv *00* outer/

# Inner domain.
python drycblles_input.py inner

# Output LBCs from parent are named `lbc_name_edge_out.time`.
# Link them here without the `_out` part.
for file in outer/lbc_*_out.0*; do
  link_name="$(basename "${file/_out/}")"
  ln -s "$PWD/$file" "$link_name"
done

mpiexec -n 8 ./microhh init drycblles
mpiexec -n 8 ./microhh run drycblles
python cross_to_nc.py -n 6
mv *.nc inner/
mv *00* inner/
