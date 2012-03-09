#ifndef DNS
#define DNS

#include "grid.h"
#include "fields.h"

class cdns
{
  public:
    cdns(cgrid *, cfields *);
    ~cdns();

    bool loop;

  private:
    cgrid *grid;
    cfields *fields;
};
#endif
