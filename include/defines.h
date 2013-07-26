#include <climits>
#define restrict __restrict__

#define kappa  0.4

#define dsmall 1.e-9
#define dbig   1.e9
#define dhuge  1.e30

#define ulhuge ULONG_MAX

// finite difference coefficients

// 4th order interpolation
#define ci0 -1./16.
#define ci1  9./16.
#define ci2  9./16.
#define ci3 -1./16.

#define bi0  5./16.
#define bi1 15./16.
#define bi2 -5./16.
#define bi3  1./16.

#define ti0  1./16.
#define ti1 -5./16.
#define ti2 15./16.
#define ti3  5./16.

// 4th order gradient
#define cg0   1.
#define cg1 -27.
#define cg2  27.
#define cg3  -1.
#define cgi   1./24.

#define bg0 -23.
#define bg1  21.
#define bg2   3.
#define bg3  -1.

#define tg0   1.
#define tg1  -3.
#define tg2 -21.
#define tg3  23.

// 4th order divgrad
#define cdg0 -1460./576.
#define cdg1   783./576.
#define cdg2   -54./576.
#define cdg3     1./576.

