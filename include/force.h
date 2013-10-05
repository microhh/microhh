#ifndef FORCE
#define FORCE

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cforce
{
  public:
    cforce(cgrid *, cfields *, cmpi *);
    ~cforce();
    int readinifile(cinput *);
    int init();
    int create(cinput *);
    int exec(double);

    std::vector<std::string> lslist;
    std::map<std::string, double *> lsprofs;

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    bool allocated;

    std::string swforce;
    std::string swls;
    std::string swwls;

    double uflow;
    double fc;

    double *ug;
    double *vg;

    double *wls;

    int flux(double * const, const double * const,
             const double * const, const double);
    int coriolis_2nd(double * const, double * const,
                     const double * const, const double * const,
                     const double * const, const double * const);
    int coriolis_4th(double * const, double * const,
                     const double * const, const double * const,
                     const double * const, const double * const);
    int lssource(double * const, const double * const);
    int advecwls_2nd(double * const, const double * const,
                     const double * const, const double * const);

    inline double interp2(const double, const double); ///< 2nd order interpolation function.
};
#endif
