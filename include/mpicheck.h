#ifndef MPICHECK
#define MPICHECK

#include "grid.h"
#include "fields.h"

class cmpicheck
{
  public:
    cmpicheck(cgrid *, cfields *);
    ~cmpicheck();
    
    int create();

  private:
    cgrid *grid;
    cfields *fields;
};
#endif
