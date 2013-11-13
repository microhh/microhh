#ifndef STATS_DNS
#define STATS_DNS

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "stats.h"

class cstats_dns : public cstats
{
  public:
    cstats_dns(cgrid *, cfields *, cmpi *);
    ~cstats_dns();

    int readinifile(cinput *);
    int init(double);
    int create(int);
    int exec(int, double, unsigned long);

  private:
    bool allocated;
    bool initialized;

    NcFile *dataFile;
    NcDim  *z_dim, *zh_dim, *t_dim;
    NcVar  *z_var, *zh_var, *t_var, *iter_var;

    double *uabs, *vabs;

    profmap profs;
    int addprof(std::string, std::string);

    int calcmean     (double *, double *, double);
    int calcmoment   (double *, double *, double *, double, int);
    int calcflux     (double *, double *, double *, double *, int, int);
    int calcdiff     (double *, double *, double *, double);
    int calcgrad     (double *, double *, double *);
    int addfluxes    (double *, double *, double *);
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
