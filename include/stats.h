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
    int init(double);
    int create(int);
    unsigned long gettimelim(unsigned long);
    int exec(int, double);
    int dostats(int, unsigned long );

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    cstats_dns *stats_dns;
    cstats_les *stats_les;

    std::string swstats;

    double statstime;
    unsigned long istatstime;

};
#endif

