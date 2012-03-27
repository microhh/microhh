#ifndef TIMELOOP
#define TIMELOOP

#include <sys/time.h>
#include "grid.h"
#include "fields.h"

class ctimeloop
{
  public:
    ctimeloop(cgrid *, cfields *);
    ~ctimeloop();

    int timestep();
    int settimestep(double);

    bool loop;
    bool adaptivestep;

    double time;
    double runtime;
    double dt;
    double cflmax;

    int iteration;
    int itime;
    int iruntime;
    int idt;


  private:
    cgrid *grid;
    cfields *fields;

    timeval start;
    timeval end;
};
#endif
