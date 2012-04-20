#ifndef FIELDS
#define FIELDS

#include "grid.h"
#include "field3d.h"

class cfields
{
  public:
    // functions
    cfields(cgrid *);
    ~cfields();

    int readinifile(cinput *);
    int initfields();
    int createfields();
    int save(int);
    int load(int);

    // int resettend();
    int boundary();

    double check(int);
    inline double interp2(const double, const double);

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
    bool allocated;

    // functions
    double calcmom (double *, double *, double *, double *);
    double calctke (double *, double *, double *, double *);
    double calcmass(double *, double *, double *, double *, double *);
};
#endif

