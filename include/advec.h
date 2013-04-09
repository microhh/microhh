#ifndef ADVEC
#define ADVEC

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

/**
 * Base class for the advection scheme.
 * This class handles the case when advection is turned off. Derived classes are
 * implemented that handle different advection schemes.
 */
class cadvec
{
  public:
    cadvec(cgrid *, cfields *, cmpi *); ///< Constructor of the advection class.
    virtual ~cadvec();                  ///< Destructor of the advection class.

    virtual int readinifile(cinput *);  ///< Processes the data from the input file.

    virtual double getcfl(double); ///< Retrieve the CFL number.
    virtual unsigned long gettimelim(unsigned long, double);
    virtual int exec();

    double cflmax;
  private:
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.
    cmpi    *mpi;    ///< Pointer to mpi class.

};
#endif
