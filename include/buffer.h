#ifndef BUFFER
#define BUFFER

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cbuffer
{
  public:
    cbuffer(cgrid *, cfields *, cmpi *);
    ~cbuffer();

    int readinifile(cinput *);
    int setbuffers();
    int init();
    int exec();
    int save();
    int load();

    // CvH make private later
    int ibuffer;

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double buffersigma;
    double bufferbeta;

    int bufferkstart;
    int bufferkcells;

    std::map<std::string, double*> bufferprofs;

    bool allocated;

    int setbuffer(double *, double *);
    int buffer   (double *, double *, double *, double *);
};
#endif
