#ifndef GRID
#define GRID

#include "input.h"

class cgrid
{
  public:
    // functions
    cgrid();
    ~cgrid();
    int readinifile(cinput *);
    int init(int, int);
    int create();
    int calculate();
    int save(int);
    int load(int);

    // variables
    int itot;
    int jtot;
    int ktot;

    int imax;
    int jmax;
    int kmax;

    int iblock;
    int kblock;

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

  private:
    bool allocated;
};
#endif
