#ifndef FIELDS
#define FIELDS
#include <map>
#include "grid.h"
#include "mpiinterface.h"
#include "field3d.h"

typedef std::map<std::string, cfield3d *> fieldmap;

class cfields
{
  public:
    // functions
    cfields(cgrid *, cmpi *);
    ~cfields();
    int readinifile(cinput *);
    int init();
    int create(cinput *);

    int initmomfld(cfield3d*&, cfield3d*&, std::string);
    int initpfld(std::string);
    int initdfld(cfield3d*&, std::string);
    
    int save(int);
    int load(int);

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
    
    fieldmap s;
    fieldmap sd;
    fieldmap sp;
    fieldmap st;
    fieldmap m;
    fieldmap mp;
    fieldmap mt;
    
    // temporary arrays
    cfield3d *tmp1;
    cfield3d *tmp2;

    // eddy viscosity for LES
    cfield3d *evisc;
    double tPr;

    double visc;
    // double viscs;
    
  private:
    // variables
    cgrid *grid;
    cmpi  *mpi;
    bool allocated;

    // perturbations
    double rndamp;
    double rndamps;
    double rndz;
    double rndbeta;
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

