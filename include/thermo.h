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

    // interfacint functions to get buoyancy properties from other classes
    virtual int getbuoyancysurf   (cfield3d *);
    virtual int getbuoyancyfluxbot(cfield3d *);
    virtual int getbuoyancy(cfield3d *, cfield3d *);

    std::string getname();

  protected:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    std::string swthermo;
};
#endif
