#ifndef BOUNDARY
#define BOUNDARY

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

struct field3dbc
{
  double bot;
  double top;
  int bcbot;
  int bctop;
};

class cboundary
{
  public:
    cboundary(cgrid *, cfields *, cmpi *);
    ~cboundary();

    int readinifile(cinput *);
    int setvalues();
    int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    int setbc      (double *, double *, double *, int, double, double);
    int setbc_patch(double *, double, double, double);

    // int setgcbot_2nd(double *, double *, int, double);
    // int setgctop_2nd(double *, double *, int, double);
    // int setgcbot_4th(double *, double *, int, double);
    // int setgctop_4th(double *, double *, int, double);

    int setgcbot_2nd(double *, double *, int, double *, double *);
    int setgctop_2nd(double *, double *, int, double *, double *);
    int setgcbot_4th(double *, double *, int, double *, double *);
    int setgctop_4th(double *, double *, int, double *, double *);

    int setgcbotw_4th(double *);
    int setgctopw_4th(double *);

    std::string swboundary;
    std::string swboundarytype;

    int mbcbot;
    int mbctop;

    typedef std::map<std::string, field3dbc *> bcmap;
    bcmap sbc;

    // patch type
    int    patch_dim;
    double patch_xh;
    double patch_xr;
    double patch_xi;
    double patch_facr;
    double patch_facl;

    inline double grad4x(const double, const double, const double, const double);
};
#endif
