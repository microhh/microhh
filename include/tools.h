#ifndef TOOLS
#define TOOLS

/* CUDA error checking, from: http://choorucode.com/2011/03/02/how-to-do-error-checking-in-cuda/
   In debug mode, CUDACHECKS is defined and all kernel calls are checked with cudaCheckError().
   All CUDA api calls are always checked with cudaSafeCall() */
#define cudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define cudaCheckError()  __cudaCheckError(__FILE__, __LINE__)

enum ReduceType {sumType, maxType}; ///< Enumerator holding the different reduction types
const int reduceMaxThreads = 0.;    ///< Maximum number of threads used in reduce algorithms

void reduceInterior(double *, double *, int, int, int, int, int, int, int, int, int, int, ReduceType);
void reduceAll(double *, double *, int, int, int, ReduceType, double);

// Wrapper to check for errors in CUDA api calls (e.g. cudaMalloc)
inline void __cudaSafeCall(cudaError err, const char *file, const int line)
{
  if (cudaSuccess != err)
  {
    printf("cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    throw 1;
  }
  return;
}
 
// Function to check for errors in CUDA kernels. Call directly after kernel.
inline void __cudaCheckError(const char *file, const int line)
{
#ifdef CUDACHECKS
  cudaError err = cudaGetLastError();
  if (cudaSuccess != err)
  {
    printf("cudaCheckError() failed at %s:%i : %s\n", file, line, cudaGetErrorString( err ) );
    throw 1;
  }
 
  err = cudaDeviceSynchronize();
  if(cudaSuccess != err)
  {
    printf("cudaCheckError() with sync failed at %s:%i : %s\n", file, line, cudaGetErrorString( err ) );
    throw 1;
  }
#endif
  return;
}

#endif
