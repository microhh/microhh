#ifndef TIMELOOP
#define TIMELOOP

#include <sys/time.h>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include <string>

class ctimeloop
{
  public:
    ctimeloop(cgrid *, cfields *, cmpi *);
    ~ctimeloop();

    int readinifile(cinput *);
    int timestep();
    int postprocstep();
    int settimestep();

    int exec();
    bool insubstep();
    double getsubdt();

    int settimelim();

    unsigned long settimelim(unsigned long);

    int save(int starttime);
    int load(int starttime);

    int docheck();
    double check();

    int dosave();

    // variables
    int substep;
    bool loop;
    bool adaptivestep;

    double time;
    double dt;
    double dtmax;
    double runtime;
    double postproctime;
    double savetime;
    double starttime;
    double dtlim;

    int istarttime;
    int iteration;

    unsigned long itime;
    unsigned long iruntime;
    unsigned long idt;
    unsigned long idtmax;
    unsigned long ipostproctime;
    unsigned long isavetime;
    unsigned long idtlim;

    double ifactor;
    double precision;

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    timeval start;
    timeval end;

    int rkorder;

    int outputiter;

    int rk3(double *, double *, double);
    int rk4(double *, double *, double);
    double rk3subdt(double);
    double rk4subdt(double);
};
#endif
