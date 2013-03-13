#ifndef TIMELOOP
#define TIMELOOP

#include <sys/time.h>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class ctimeloop
{
  public:
    ctimeloop(cgrid *, cfields *, cmpi *);
    ~ctimeloop();

    int readinifile(cinput *);
    int timestep();
    int postprocstep();
    int settimestep(double, double);

    int exec();
    bool insubstep();
    double getsubdt();

    int save(int);
    int load(int);

    int docheck();
    double check();

    int dosave();
    int dostats();

    // variables
    int substep;
    bool loop;
    bool adaptivestep;

    double time;
    double dt;
    double dtmax;
    double cflmax;
    double dnmax;
    double runtime;

    int iteration;

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double ifactor;

    timeval start;
    timeval end;

    int rkorder;

    unsigned long itime;
    unsigned long iruntime;
    unsigned long idt;
    unsigned long idtmax;

    int startiter;
    int outputiter;
    int saveiter;
    int statsiter;
    int maxiter;
    int postprociter;

    int rk3(double *, double *, double);
    int rk4(double *, double *, double);
    double rk3subdt(double);
    double rk4subdt(double);
};
#endif
