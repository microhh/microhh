#ifndef FIELDS
#define FIELDS

#include "grid.h"
#include "mpiinterface.h"
#include "field3d.h"

class cfields
{
  public:
    // functions
    cfields(cgrid *, cmpi *);
    ~cfields();

    int readinifile(cinput *);
    int init();
    int create(cinput *);
    int save(int);
    int load(int);

    // int resettend();
    int boundary();

    double checkmom ();
    double checktke ();
    double checkmass();

    // variables
    cfield3d *u;
    cfield3d *v;
    cfield3d *w;
    cfield3d *p;

    cfield3d *ut;
    cfield3d *vt;
    cfield3d *wt;

    cfield3d *s;
    cfield3d *st;

    double visc;
    double viscs;

  private:
    // variables
    cgrid *grid;
    cmpi  *mpi;
    bool allocated;

    // perturbations
    double rndamp;
    double vortexamp;
    int nvortexpair;
    int vortexaxis;

    // functions
    double calcmom_2nd(double *, double *, double *, double *);
    double calctke_2nd(double *, double *, double *, double *);
    double calcmass   (double *, double *);
    inline double interp2(const double, const double);
};
#endif

