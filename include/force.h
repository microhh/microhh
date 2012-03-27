#ifndef FORCE
#define FORCE

#include "grid.h"
#include "fields.h"

class cforce
{
  public:
    cforce(cgrid *, cfields *);
    ~cforce();

    int exec(double);

  private:
    cgrid *grid;
    cfields *fields;
    int flux(double *, double *, double *, double);

};
#endif
