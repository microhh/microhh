#ifndef GRID
#define GRID
class grid
{
  public:
    // functions
    grid();
    ~grid();

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
    
    double *x;
    double *y;
    double *z;
    double *xh;
    double *yh;
    double *zh;
};
#endif
