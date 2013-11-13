#ifndef STATS
#define STATS

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
// #include "stats_dns.h"
// #include "stats_les.h"

struct statsvar
{
  NcVar *ncvar;
  double *data;
};
typedef std::map<std::string, statsvar> profmap;

class cstats
{
  public:
    cstats(cgrid *, cfields *, cmpi *);
    virtual ~cstats();

    virtual int readinifile(cinput *);
    virtual int init(double);
    virtual int create(int);
    virtual int exec(int, double, unsigned long);

    unsigned long gettimelim(unsigned long);
    int dostats(int, unsigned long);

  protected:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double statstime;
    unsigned long istatstime;
};
#endif

