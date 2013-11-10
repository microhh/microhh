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
    int initdfld(std::string);
    
    int save(int);
    int load(int);

    double checkmom ();
    double checktke ();
    double checkmass();

    // 3d fields for momentum
    cfield3d *u;
    cfield3d *v;
    cfield3d *w;

    cfield3d *ut;
    cfield3d *vt;
    cfield3d *wt;

    // maps of 3d fields
    fieldmap ap;
    fieldmap at;

    fieldmap mp;
    fieldmap mt;

    fieldmap s;
    fieldmap sd;
    fieldmap sp;
    fieldmap st;

    // TODO remove these to and bring them to diffusion model
    double tPr;
    double visc;

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
    int addmeanprofile(cinput *, std::string, double *, double);
    int randomnize(cinput *, std::string, double *);
    int addvortexpair(cinput* inputin);
    double calcmass(double *, double *);
    inline double interp2(const double, const double);
};
#endif

