Simple example case of how to setup the Immersed Boundaries for
objects defined as polygons.

The polygons are defined in `poly_coordinates.txt`. To run the case,
first create the vertical grid and input profiles:

    python ib_poly_prof.py
    ./microhh init ib_poly

Now that the 3D grid is known (file `grid.0000000`), the IB can be
initialized:

    python ib_poly_init.py

This creates several text files (`*.ib_input`) which define the
IB properties (location ghost cells, boundary, interpolation points, ...).

Finally, the case can be started with:

    ./microhh run ib_poly
