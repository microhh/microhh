#ifndef MODEL
#define MODEL

#include <string>

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "boundary.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "buoyancy.h"
#include "pres.h"
#include "buffer.h"
#include "timeloop.h"
#include "stats.h"
#include "cross.h"

class cmodel
{
  public:
    cmodel(cgrid *, cmpi *);
    ~cmodel();
    int readinifile(cinput *);
    int init();
    int load();
    int save(cinput *);
    int exec();

  private:
    cgrid   *grid;
    cmpi    *mpi;

    // switches for included schemes
    int iadvec;

    // fields to be created
    cfields   *fields;

    // model operators
    cboundary *boundary;
    ctimeloop *timeloop;
    cadvec    *advec;
    cdiff     *diff;
    cpres     *pres;  
    cforce    *force;   
    cbuoyancy *buoyancy;
    cbuffer   *buffer;

    // load the postprocessing modules
    cstats    *stats;
    ccross    *cross;
};
#endif
