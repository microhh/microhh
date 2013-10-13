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

    std::vector<std::string> o1;
    std::vector<std::string> o2;
    std::vector<std::string> o3;
    std::vector<std::string> o4;
    std::vector<std::string> flux;
    std::vector<std::string> diff;

    struct statsvar
    {
      std::string name;
      NcVar *ncvar;
      double *data;
    };

    std::map<std::string, statsvar> statsvarmap;

    NcFile *dataFile;
    NcDim  *z_dim, *zh_dim, *t_dim;
    NcVar  *z_var, *zh_var, *t_var, *iter_var;
    statsvar  u, v, w, s;
    statsvar  evisc;
    statsvar  u2, v2, w2, s2;
    statsvar  u3, v3, w3, s3;
    statsvar  ugrad, vgrad, sgrad;
    statsvar  wu, wv, ws;
    statsvar  udiff, vdiff, sdiff;
    statsvar  uflux, vflux, sflux;

    int calcmean     (double *, double *);
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
