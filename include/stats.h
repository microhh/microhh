#ifndef STATS
#define STATS

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cstats
{
  public:
    cstats(cgrid *, cfields *, cmpi *);
    ~cstats();

    int readinifile(cinput *);
    int init();
    int exec(int, double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int istats;

    NcFile *dataFile;
    NcDim  *z_dim, *zh_dim, *t_dim;
    NcVar  *z_var, *zh_var, *t_var;
    NcVar  *u_var , *v_var , *w_var , *s_var ;
    NcVar  *u2_var, *v2_var, *w2_var, *s2_var;
    NcVar  *wu_var, *wv_var, *ws_var;
    NcVar  *wud_var, *wvd_var, *wsd_var;

    double *u , *v , *w , *s ;
    double *u2, *v2, *w2, *s2;
    double *wu , *wv , *ws ;
    double *wud, *wvd, *wsd;

    int calcmean(double *, double *, int);
    int calcvar (double *, double *, double *, int);
    int calcflux(double *, double *, double *);
    int calcdiff(double *, double *, double *, double);
    
    int nstats;
};
#endif

