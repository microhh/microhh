# Ubuntu 14.04
if(USEMPI) 
  set(ENV{CC}  mpicc ) # C compiler for parallel build
  set(ENV{CXX} mpicxx) # C++ compiler for parallel build
else()
  set(ENV{CC}  icc ) # C compiler for serial build
  set(ENV{CXX} icpc) # C++ compiler for serial build
endif()

set(USER_CXX_FLAGS "")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib/x86_64-linux-gnu/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_C       "/usr/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/usr/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/usr/lib/x86_64-linux-gnu/libhdf5.a")
set(HDF5_LIB_2         "/usr/lib/x86_64-linux-gnu/libhdf5_hl.a")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=__restrict__)
