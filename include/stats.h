#ifndef STATS
#define STATS

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "stats_dns.h"
#include "stats_les.h"

class cstats
{
  public:
    cstats(cgrid *, cfields *, cmpi *);
    ~cstats();

    int readinifile(cinput *);
    int init();
    int create(int);
    int exec(int, double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    cstats_dns *stats_dns;
    cstats_les *stats_les;

    int istats;
};
#endif

