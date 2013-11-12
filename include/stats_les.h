#ifndef STATS_LES
#define STATS_LES

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cstats_les
{
  public:
    cstats_les(cgrid *, cfields *, cmpi *);
    ~cstats_les();

    int readinifile(cinput *);
    int init();
    int create(int);
    int exec(int, double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    bool allocated;
    bool initialized;

    NcFile *dataFile;
    NcDim  *z_dim, *zh_dim, *t_dim;
    NcVar  *z_var, *zh_var, *t_var, *iter_var;
    NcVar  *u_var , *v_var, *w_var;
    NcVar  *evisc_var;
    NcVar  *u2_var, *v2_var, *w2_var;
    NcVar  *u3_var, *v3_var, *w3_var, *s3_var;
    NcVar  *ugrad_var, *vgrad_var, *sgrad_var;
    NcVar  *wu_var, *wv_var, *ws_var;
    NcVar  *udiff_var, *vdiff_var, *sdiff_var;
    NcVar  *uflux_var, *vflux_var, *sflux_var;

    double *u , *v , *w;
    double *uabs, *vabs;
    double *evisc;
    double *u2, *v2, *w2;
    double *u3, *v3, *w3, *s3;
    double *wu , *wv , *ws;
    double *ugrad, *vgrad, *sgrad;
    double *udiff, *vdiff, *sdiff;
    double *uflux, *vflux, *sflux;

    struct statsvar
    {
      NcVar *ncvar;
      double *data;
    };

    std::map<std::string, statsvar> profs;
    int addprof(std::string, std::string);

    int calcmean     (double *, double *, double);
    int calcmoment   (double *, double *, double *, double, int);
    int calcdiff     (double *, double *, double *, double *, double *, double *, double);
    int calcgrad     (double *, double *, double *);
    int calcflux     (double *, double *, double *, double *, int, int);
    int calctkebudget(double *, double *, double *, double *, double *,
                      double *, double *,
                      double *, double *,
                      double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *,
                      double *, double *, double *,
                      double *, double *,
                      double *, double *, double);

    int nstats;
};
#endif
