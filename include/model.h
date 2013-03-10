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
    cmodel(cgrid *, cfields *, cmpi *, std::string);
    ~cmodel();
    int readinifile(cinput *);
    int init();
    int load();
    int save();
    int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
    std::string simname;

    // create the boundary conditions class
    cboundary *boundary;

    // create the instances of the model operations
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
