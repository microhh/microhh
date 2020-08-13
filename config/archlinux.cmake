# ARCH Linux / Manjaro
if(USEMPI) 
  set(ENV{CC}  mpicc)    # C compiler for parallel build
  set(ENV{CXX} mpicxx)   # C++ compiler for parallel build
  set(ENV{FC}  mpif90)   # Fortran compiler for parallel build
else()
  # NOTE: CUDA is picky on which GCC version is used.
  # See: https://gist.github.com/ax3l/9489132#user-content-nvcc
  # GCC <= 8 is required for CUDA 10.2.xx
  set(ENV{CC}  gcc-8)      # C compiler for serial build
  set(ENV{CXX} g++-8)      # C++ compiler for serial build
  set(ENV{FC}  gfortran-8) # Fortran compiler for serial build
endif()

set(USER_CXX_FLAGS "-std=c++14 -fopenmp")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib/libfftw3.so")
set(FFTWF_LIB          "/usr/lib/libfftw3f.so")
set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_C       "/usr/lib/libnetcdf.so")
set(NETCDF_LIB_CPP     "/usr/lib/libnetcdf_c++4.so")
set(HDF5_LIB_1         "/usr/lib/libhdf5.a")
set(HDF5_LIB_2         "/usr/lib/libhdf5_hl.a")
set(SZIP_LIB           "/usr/lib/libsz.so")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl dl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(CUFFT_LIB "/opt/cuda/lib64/libcufft.so")
  set(LIBS ${LIBS} ${CUFFT_LIB} -rdynamic )
  set(USER_CUDA_NVCC_FLAGS "-arch=sm_70")
  list(APPEND CUDA_NVCC_FLAGS " -std=c++14")
  set(USER_CUDA_NVCC_FLAGS_RELEASE "-Xptxas -O3 -use_fast_math")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
