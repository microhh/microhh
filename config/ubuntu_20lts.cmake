# Ubuntu 20.04
if(USEMPI) 
  set(ENV{CC}  mpicc ) # C compiler for parallel build
  set(ENV{CXX} mpicxx) # C++ compiler for parallel build
else()
  set(ENV{CC}  gcc) # C compiler for serial build
  set(ENV{CXX} g++) # C++ compiler for serial build
endif()

if(USECUDA)
  set(USER_CXX_FLAGS "-std=c++14 -fopenmp")
else()
  set(USER_CXX_FLAGS "-std=c++14")
endif()

set(USER_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")
set(USER_FC_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib/x86_64-linux-gnu/libfftw3.so")
set(FFTWF_LIB          "/usr/lib/x86_64-linux-gnu/libfftw3f.so")
set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_C       "/usr/lib/x86_64-linux-gnu/libnetcdf.so")
set(HDF5_LIB_1         "/usr/lib/x86_64-linux-gnu/libhdf5_serial.so")
set(HDF5_LIB_2         "/usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
  set(USER_CUDA_NVCC_FLAGS "-std=c++14 -arch=sm_70")
  set(USER_CUDA_NVCC_FLAGS_RELEASE "-Xptxas -O3")
  set(USER_CUDA_FLAGS_DEBUG "-Xptxas -O3 -DCUDACHECKS")
  set(LIBS ${LIBS} cufft)
  add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_OWN)
  #add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_CUDA)
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
add_definitions(-DRTE_RRTMGP_USE_CBOOL)
