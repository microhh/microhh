#ifndef TOOLS
#define TOOLS

/* CUDA error checking, from: http://choorucode.com/2011/03/02/how-to-do-error-checking-in-cuda/
   In debug mode, CUDACHECKS is defined and all calls are checked. If not, the calls shouldn't
   have any overhead. */
#define cudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define cudaCheckError()  __cudaCheckError(__FILE__, __LINE__)

void reduceInterior(double *, double *, int, int, int, int, int, int, int, int, int, int, int);
void reduceAll(double *, double *, int, int, int, int, double);

// Wrapper to check for errors in CUDA api calls (e.g. cudaMalloc)
inline void __cudaSafeCall(cudaError err, const char *file, const int line)
{
#ifdef CUDACHECKS
  if (cudaSuccess != err)
  {
    printf("cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    throw 1;
  }
#endif
  return;
}
 
// Wrapper to check for errors in CUDA kernels
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
