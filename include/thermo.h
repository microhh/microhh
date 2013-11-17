#ifndef THERMO
#define THERMO

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cthermo
{
  public:
    cthermo(cgrid *, cfields *, cmpi *);
    virtual ~cthermo();
    virtual int readinifile(cinput *);
    virtual int exec();

  protected:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double gravitybeta; // gravity multiplied with thermal expansion coefficient
};
#endif
