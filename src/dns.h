#ifndef DNS
#define DNS

#include "grid.h"
#include "fields.h"

class cdns
{
  public:
    cdns(cgrid *, cfields *);
    ~cdns();

    int timestep();

    bool loop;
    double time;
    double runtime;
    double dt;

    int iteration;
    int itime;
    int iruntime;
    int idt;


  private:
    cgrid *grid;
    cfields *fields;
};
#endif
