#ifndef TIMELOOP
#define TIMELOOP

#include <sys/time.h>
#include "input.h"
#include "grid.h"
#include "fields.h"

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

    int save(int);
    int load(int);

    // variables
    int substep;
    bool loop;
    bool adaptivestep;

    double time;
    double dt;
    double cflmax;
    double runtime;

    int iteration;
    int itime;
    int iruntime;
    int idt;
    int maxiter;

  private:
    cgrid *grid;
    cfields *fields;

    double ifactor;

    timeval start;
    timeval end;

    int rkorder;

    int rk3(double *, double *, double);
    int rk4(double *, double *, double);
    double rk3subdt(double);
    double rk4subdt(double);
};
#endif
