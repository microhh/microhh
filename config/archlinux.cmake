# ARCH Linux
#set(ENV{CXX} g++) # compiler for serial build
set(ENV{CC}  mpicc)
set(ENV{CXX} mpicxx) # compiler for parallel build

set(USER_CXX_FLAGS "")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib64/libfftw3.so")
set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_C       "/usr/lib64/libnetcdf.so")
set(NETCDF_LIB_CPP     "/usr/lib64/libnetcdf_c++.so")
set(HDF5_LIB_1         "/usr/lib64/libhdf5.a")
set(HDF5_LIB_2         "/usr/lib64/libhdf5_hl.a")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(LIBS ${LIBS} -rdynamic /opt/cuda/lib64/libcufft.so)
  set(USER_CUDA_NVCC_FLAGS "-arch=sm_21")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
