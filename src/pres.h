#ifndef PRES
#define PRES

#include "grid.h"
#include "fields.h"

class cpres
{
  public:
    cpres(cgrid *, cfields *);
    ~cpres();

    int exec();

  private:
    cgrid *grid;
    cfields *fields;
};
#endif
