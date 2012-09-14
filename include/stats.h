#ifndef STATS
#define STATS

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cstats
{
  public:
    cstats(cgrid *, cfields *, cmpi *);
    ~cstats();
    int readinifile(cinput *);
    int exec(int);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int istats;
};
#endif

