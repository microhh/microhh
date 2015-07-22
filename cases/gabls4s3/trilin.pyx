#cython: boundscheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as np
DTYPEd = np.double
ctypedef np.double_t DTYPEd_t
DTYPEi = np.int
ctypedef np.int_t DTYPEi_t

def trilin_python(field_in, field_out, xi, yi, zi, xo, yo, zo):
    for i in range(xi.size):
        print(i, xi.size)
        ii = xi[i]  
        xb = xo[i]
        xa = 1.-xb
        for j in range(yo.size):
            jj = yi[j]
            yb = yo[j]
            ya = 1.-yb
            for k in range(zo.size):
                kk = zi[k]
                zb = zo[k]
                za = 1.-zb
     
                field_out[k,j,i] = za * ya * xa * field_in[kk,   jj,   ii  ] + \
                                   zb * ya * xa * field_in[kk+1, jj,   ii  ] + \
                                   za * yb * xa * field_in[kk,   jj+1, ii  ] + \
                                   za * ya * xb * field_in[kk,   jj,   ii+1] + \
                                   zb * yb * xa * field_in[kk+1, jj+1, ii  ] + \
                                   zb * ya * xb * field_in[kk+1, jj,   ii+1] + \
                                   za * yb * xb * field_in[kk,   jj+1, ii+1] + \
                                   zb * yb * xb * field_in[kk+1, jj+1, ii+1]

def trilin_cython(np.ndarray[DTYPEd_t, ndim=3] field_in, \
                  np.ndarray[DTYPEd_t, ndim=3] field_out, \
                  np.ndarray[DTYPEi_t, ndim=1] xi, \
                  np.ndarray[DTYPEi_t, ndim=1] yi, \
                  np.ndarray[DTYPEi_t, ndim=1] zi, \
                  np.ndarray[DTYPEd_t, ndim=1] xo, \
                  np.ndarray[DTYPEd_t, ndim=1] yo, \
                  np.ndarray[DTYPEd_t, ndim=1] zo):

    cdef int i,j,k,ii,jj,kk
    cdef double xa,xb,ya,yb,za,zb

    cdef int nx = xi.size
    cdef int ny = yi.size
    cdef int nz = zi.size

    for i in range(nx):
        ii = xi[i]  
        xb = xo[i]
        xa = 1.-xb
        for j in range(ny):
            jj = yi[j]
            yb = yo[j]
            ya = 1.-yb
            for k in range(nz):
                kk = zi[k]
                zb = zo[k]
                za = 1.-zb
     
                field_out[k,j,i] = za * ya * xa * field_in[kk,   jj,   ii  ] + \
                                   zb * ya * xa * field_in[kk+1, jj,   ii  ] + \
                                   za * yb * xa * field_in[kk,   jj+1, ii  ] + \
                                   za * ya * xb * field_in[kk,   jj,   ii+1] + \
                                   zb * yb * xa * field_in[kk+1, jj+1, ii  ] + \
                                   zb * ya * xb * field_in[kk+1, jj,   ii+1] + \
                                   za * yb * xb * field_in[kk,   jj+1, ii+1] + \
                                   zb * yb * xb * field_in[kk+1, jj+1, ii+1]
