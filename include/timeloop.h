#ifndef TIMELOOP
#define TIMELOOP

#include <sys/time.h>
#include "input.h"
#include "grid.h"
#include "fields.h"

#define restrict __restrict__

class ctimeloop
{
  public:
    ctimeloop(cgrid *, cfields *);
    ~ctimeloop();

    int readinifile(cinput *);
    int timestep();
    int settimestep(double);

    int exec();
    bool insubstep();
    double getsubdt();

    // variables
    int substep;
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

    int rk3(double *, double *, double);
    int rk4(double *, double *, double);
    double rk3subdt(double);
    double rk4subdt(double);
};
#endif
