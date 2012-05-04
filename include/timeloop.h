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
    int settimestep(double, double);

    int exec();
    bool insubstep();
    double getsubdt();

    int save(int, int);
    int load(int, int);

    int docheck();
    double check();

    int dosave();

    // variables
    int substep;
    bool loop;
    bool adaptivestep;

    double time;
    double dt;
    double cflmax;
    double dnmax;
    double runtime;

    int iteration;

  private:
    cgrid *grid;
    cfields *fields;

    double ifactor;

    timeval start;
    timeval end;

    int rkorder;

    unsigned long itime;
    unsigned long iruntime;
    unsigned long idt;

    int outputiter;
    int saveiter;
    int maxiter;

    int rk3(double *, double *, double);
    int rk4(double *, double *, double);
    double rk3subdt(double);
    double rk4subdt(double);
};
#endif
