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
    int create(cinput *);
    int init();
    int exec();

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double bufferz;     ///< Height above which the buffer is applied
    double buffersigma;
    double bufferbeta;

    int bufferkstart;
    int bufferkstarth;
    int bufferkcells;

    std::map<std::string, double*> bufferprofs;

    bool allocated;
    std::string swbuffer;

    // int setbuffer(double *, double *);
    int buffer   (double *, double *, double *, double *);
};
#endif
