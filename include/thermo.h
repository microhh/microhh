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
    virtual int create();
    virtual int getbuoyancysurf(cfield3d *);

    std::string getname();

  protected:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    std::string swthermo;
};
#endif
