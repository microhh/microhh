#ifndef BUFFER
#define BUFFER

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cbuffer
{
  public:
    cbuffer(cgrid *, cfields *, cmpi *); ///< Constructor of the buffer class.
    ~cbuffer();                          ///< Destructor of the buffer class.

    int readinifile(cinput *); ///< Processing data of the input file.
    int init();                ///< Initialize the arrays that contain the profiles.
    int create(cinput *);      ///< Read the profiles of the forces from the input.
    int exec();                ///< Add the tendencies created by the damping.

  private:
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.
    cmpi    *mpi;    ///< Pointer to mpi class.

    double bufferz;     ///< Height above which the buffer is applied.
    double buffersigma; ///< Damping frequency.
    double bufferbeta;  ///< Exponent for damping increase with height.

    int bufferkstart;  ///< Grid point at cell center at which damping starts.
    int bufferkstarth; ///< Grid point at cell face at which damping starts.

    std::map<std::string, double*> bufferprofs; ///< Map containing the buffer profiles.

    bool allocated; ///< Boolean flag to indicate allocation of arrays.

    std::string swbuffer; ///< Switch for buffer.

    int buffer(double * const, const double * const, 
              const double * const, const double * const); ///< Calculate the tendency 
};
#endif
