Simple example case of how to setup the Immersed Boundaries for
objects defined as polygons. Requires the `microhh_tools.py` and `ib_tools.py`
scripts from the `microhh_root/python/` and `microhh_root/python/IB/` directories.

The polygons are defined in `poly_coordinates.txt`, where each line defines an 
object in the format `z_object, x1, y1, x2, y2, .., .., xn, yn`. To run the case,
first create the input profiles:

    python ib_poly_prof.py

Now that the vertical grid is defined, the IB can be initialized:

    python ib_poly_init.py

This creates several text files (`*.ib_input`) which define the
IB properties (location ghost cells, boundary, interpolation points, ...).

Finally, the case can be started with:

    ./microhh init ib_poly
    ./microhh run ib_poly

Note that when changing any of the grid properties, all steps above
have to be repeated.
