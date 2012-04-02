#ifndef GRID
#define GRID

#include "input.h"

class cgrid
{
  public:
    // functions
    cgrid(cinput *);
    ~cgrid();
    int initgrid();
    int creategrid();
    int calculate();
    int save();
    int load();

    // variables
    int itot;
    int jtot;
    int ktot;

    int imax;
    int jmax;
    int kmax;

    int igc;
    int jgc;
    int kgc;

    int icells;
    int jcells;
    int kcells;
    int ncells;
    int istart;
    int jstart;
    int kstart;
    int iend;
    int jend;
    int kend;

    double xsize;
    double ysize;
    double zsize;

    double dx;
    double dy;
    double *dz;
    double *dzh;
    double *dzi;
    double *dzhi;
    
    double *x;
    double *y;
    double *z;
    double *xh;
    double *yh;
    double *zh;
};
#endif
