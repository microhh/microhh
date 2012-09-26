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
};
#endif

