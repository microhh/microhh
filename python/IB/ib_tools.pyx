import cython
import numpy as np
cimport numpy as np
from libcpp cimport bool

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double is_left(double xp, double yp, double x0, double y0, double x1, double y1):
    return (x1-x0) * (yp-y0) - (xp-x0) * (y1-y0)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef bool is_inside(double xp, double yp,\
                    np.ndarray[np.double_t, ndim=1] x_set,\
                    np.ndarray[np.double_t, ndim=1] y_set,
                    int n):

    cdef int i,j
    cdef double wn = 0

    for i in range(n-1):
        # First check: see if point exactly on line segment:
        if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) == 0):
            return True 

        if (y_set[i] <= yp):
            if (y_set[i+1] > yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) > 0):
                    wn += 1
        else:
            if (y_set[i+1] <= yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) < 0):
                    wn -= 1

    if (wn == 0):
        return False
    else:
        return True


@cython.boundscheck(False)
@cython.wraparound(False)
cdef flag_inside(np.ndarray[np.double_t, ndim=1] xgr,\
                 np.ndarray[np.double_t, ndim=1] ygr,\
                 np.ndarray[np.double_t, ndim=1] xobj,\
                 np.ndarray[np.double_t, ndim=1] yobj,
                 np.ndarray[np.int_t,    ndim=2] inside,
                 int itot, int jtot, int idobj, int objsize):

    cdef int i,j

    for i in range(itot):
        for j in range(jtot):
            if (is_inside(xgr[i], ygr[j], xobj, yobj, objsize)):
                inside[j,i] = idobj

@cython.boundscheck(False)
@cython.wraparound(False)
def flag_inside_ib_2d(np.ndarray[np.double_t, ndim=1] x,\
                      np.ndarray[np.double_t, ndim=1] y,
                      poly):

    cdef int itot = x.size
    cdef int jtot = y.size
    cdef int n, idobj

    cdef np.ndarray[np.int_t, ndim=2] inside = np.ones([jtot, itot], dtype=np.int) * -1

    # Loop over objects inside dict outside the fast Cython code
    for idobj in range(len(poly)):
        print('Working on object {0:4d} of {1:4d}'.format(idobj+1, len(poly)))

        n = poly[idobj]['x'].size
        flag_inside(x, y, poly[idobj]['x'], poly[idobj]['y'], inside, itot, jtot, idobj, n)

    return inside
