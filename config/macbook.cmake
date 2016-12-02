# MacBook
if(USEMPI) 
  set(ENV{CC}  mpicc ) # C compiler for parallel build
  set(ENV{CXX} mpicxx) # C++ compiler for parallel build
else()
  set(ENV{CC}  gcc-mp-4.9) # C compiler for serial build
  set(ENV{CXX} g++-mp-4.9) # C++ compiler for serial build
endif()

set(GNU_SED "gsed")

set(USER_CXX_FLAGS "")
set(USER_CXX_FLAGS_RELEASE "-Ofast -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/opt/local/include")
set(FFTW_LIB           "/opt/local/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/opt/local/include")
set(NETCDF_LIB_C       "/opt/local/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/opt/local/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/opt/local/lib/libhdf5.a")
set(HDF5_LIB_2         "/opt/local/lib/libhdf5_hl.a")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(LIBS ${LIBS} -rdynamic /usr/local/cuda/lib/libcufft.dylib)
  set(USER_CUDA_NVCC_FLAGS "-arch=sm_20")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
