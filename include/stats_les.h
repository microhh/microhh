#ifndef STATS_LES
#define STATS_LES

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "stats.h"

class cstats_les : public cstats
{
  public:
    cstats_les(cgrid *, cfields *, cmpi *);
    ~cstats_les();

    int readinifile(cinput *);
    int init(double);
    int create(int);
    unsigned long gettimelim(unsigned long);
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
    int calcdiff     (double *, double *, double *, double *, double *, double *, double);
    int calcgrad     (double *, double *, double *);
    int calcflux     (double *, double *, double *, double *, int, int);
    int addfluxes    (double *, double *, double *);

    int nstats;
};
#endif
