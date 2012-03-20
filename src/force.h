#ifndef FORCE
#define FORCE

#include "grid.h"
#include "fields.h"

class cforce
{
  public:
    cforce(cgrid *, cfields *);
    ~cforce();

    int exec();

  private:
    cgrid *grid;
    cfields *fields;

};
#endif
