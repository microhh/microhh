#ifndef DIFF
#define DIFF

#include "grid.h"
#include "fields.h"

class cdiff
{
  public:
    cdiff(cgrid *, cfields *);
    ~cdiff();

    int init();
    int exec();

  private:
    cgrid *grid;
    cfields *fields;

    int diffc_2nd(double *, double *, double *, double *, double);
    int diffw_2nd(double *, double *, double *, double *, double);
};
#endif
