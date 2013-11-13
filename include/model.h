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
    cmodel(cmpi *, cinput *);
    ~cmodel();
    int readinifile();
    int init();
    int create();
    int load();
    int save();
    int exec();

  private:
    cmpi    *mpi;
    cinput  *input;

    // switches for included schemes
    std::string swadvec;
    std::string swdiff;
    std::string swpres;
    std::string swboundary;

    std::string swstats;

    cgrid   *grid;
    cfields *fields;

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

    int outputfile(bool);
    int settimestep();
};
#endif
