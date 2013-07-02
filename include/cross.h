#ifndef CROSS
#define CROSS

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class ccross
{
  public:
    ccross(cgrid *, cfields *, cmpi *);
    ~ccross();

    int readinifile(cinput *);
    int init(int);
    unsigned long gettimelim(unsigned long);
    int exec(double, int, int);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double crosstime;
    unsigned long icrosstime;

    std::string swcross;

    std::vector<int> jxz;
    std::vector<int> kxy;

    std::vector<std::string> simple;
    std::vector<std::string> lngrad;

    int crosssimple(double *, double *, std::string, int, int);
    int crosslngrad(double *, double *, double *, double *, std::string, int, int);
};
#endif

