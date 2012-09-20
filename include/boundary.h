#ifndef BOUNDARY
#define BOUNDARY

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cboundary
{
  public:
    cboundary(cgrid *, cfields *, cmpi *);
    ~cboundary();

    int readinifile(cinput *);
    int exec();

    // CvH make private later
    int iboundary;

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int setgcbot_2nd(double *, double *, int, double);
    int setgctop_2nd(double *, double *, int, double);
    int setgcbot_4th(double *, double *, int, double);
    int setgctop_4th(double *, double *, int, double);

    int setgcbotw_4th(double *);
    int setgctopw_4th(double *);

    // int iboundary;

    int bcbotmom;
    int bctopmom;

    int bcbotscal;
    int bctopscal;

    double sbot;
    double stop;

    inline double grad4x(const double, const double, const double, const double);
};
#endif
