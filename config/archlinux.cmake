# ARCH Linux / Manjaro
if(USEMPI)
  set(ENV{CC}  mpicc)    # C compiler for parallel build
  set(ENV{CXX} mpicxx)   # C++ compiler for parallel build
  set(ENV{FC}  mpif90)   # Fortran compiler for parallel build
else()
  if(USECUDA)
    # NOTE: CUDA is picky on which GCC version is used.
    # See: https://gist.github.com/ax3l/9489132#user-content-nvcc
    # GCC <= 8 is required for CUDA 10.2.xx
    set(ENV{CC}  gcc-11)      # C compiler for serial build
    set(ENV{CXX} g++-11)      # C++ compiler for serial build
    set(ENV{FC}  gfortran-11) # Fortran compiler for serial build
  else()
    set(ENV{CC}  gcc)        # C compiler for serial build
    set(ENV{CXX} g++)        # C++ compiler for serial build
    set(ENV{FC}  gfortran)   # Fortran compiler for serial build
  endif()
endif()

if(USECUDA)
  set(USER_CXX_FLAGS "-std=c++17 -fopenmp")
else()
  set(USER_CXX_FLAGS "-std=c++17")
endif()

set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
set(USER_FC_FLAGS_RELEASE "-DNDEBUG -Ofast -march=native")
set(USER_FC_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib/libfftw3.so")
set(FFTWF_LIB          "/usr/lib/libfftw3f.so")
set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_C       "/usr/lib/libnetcdf.so")
set(HDF5_LIB_1         "/usr/lib/libhdf5.so")
set(HDF5_LIB_2         "/usr/lib/libhdf5_hl.so")
set(SZIP_LIB           "/usr/lib/libsz.so")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl dl rt)

if(USECUDA)
  set(CUDA_INCLUDE "/opt/cuda/include")
  set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR} ${CUDA_INCLUDE})
else()
  set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})
endif()

#
# NEW versions
#
if(USECUDA)
  set(USER_CUDA_NVCC_FLAGS "-std=c++17 -arch=sm_70")
  set(USER_CUDA_NVCC_FLAGS_RELEASE "-Xptxas -O3")
  set(USER_CUDA_NVCC_FLAGS_DEBUG "-Xptxas -O0 -g -DCUDACHECKS")
  set(LIBS ${LIBS} cufft)
  add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_OWN)
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
add_definitions(-DRTE_RRTMGP_USE_CBOOL)

#
# OLD VERSIONS:
#
#if(USECUDA)
#  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
#  set(CUFFT_LIB "/opt/cuda/lib64/libcufft.so")
#  set(LIBS ${LIBS} ${CUFFT_LIB} -rdynamic)
#  set(USER_CUDA_NVCC_FLAGS "-arch=sm_70 -std=c++14 -Xcompiler -fopenmp")
#  set(USER_CUDA_NVCC_FLAGS_RELEASE "-Xptxas -O3")
#  set(USER_CUDA_NVCC_FLAGS_DEBUG "-O0 -g -G")
#endif()
#
#add_definitions(-DRESTRICTKEYWORD=__restrict__)
#add_definitions(-DUSE_CBOOL)
