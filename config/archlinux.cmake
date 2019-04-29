# ARCH Linux / Manjaro
if(USEMPI) 
  set(ENV{CC}  mpicc)    # C compiler for parallel build
  set(ENV{CXX} mpicxx)   # C++ compiler for parallel build
  set(ENV{FC}  mpif90)   # Fortran compiler for parallel build
else()
  set(ENV{CC}  gcc)      # C compiler for serial build
  set(ENV{CXX} g++)      # C++ compiler for serial build
  set(ENV{FC}  gfortran) # Fortran compiler for serial build
endif()

set(USER_CXX_FLAGS "-std=c++14")
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
  set(LIBS ${LIBS} -rdynamic /opt/cuda/lib64/libcufft.so)
  set(USER_CUDA_NVCC_FLAGS "-arch=sm_21")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
