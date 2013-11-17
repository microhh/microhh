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
    virtual ~cboundary();

    virtual int readinifile(cinput *);
    virtual int init();
    virtual int setvalues();

    virtual int save(int);
    virtual int load(int);

    int exec();

  protected:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    std::string swspatialorder;

    int mbcbot;
    int mbctop;

    typedef std::map<std::string, field3dbc *> bcmap;
    bcmap sbc;

    int setbc(double *, double *, double *, int, double, double, double);
    
  private:
    virtual int bcvalues();

    int setgcbot_2nd(double *, double *, int, double *, double *);
    int setgctop_2nd(double *, double *, int, double *, double *);
    int setgcbot_4th(double *, double *, int, double *, double *);
    int setgctop_4th(double *, double *, int, double *, double *);

    int setgcbotw_4th(double *);
    int setgctopw_4th(double *);

    inline double grad4x(const double, const double, const double, const double);
};
#endif
