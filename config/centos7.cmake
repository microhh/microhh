# CentOS 7.x
if(USEMPI) 
  set(ENV{CC}  mpicc ) # C compiler for parallel build
  set(ENV{CXX} mpicxx) # C++ compiler for parallel build
else()
  set(ENV{CC}  gcc) # C compiler for serial build
  set(ENV{CXX} g++) # C++ compiler for serial build
endif()

set(USER_CXX_FLAGS "-std=c++11")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR      "/usr/include")
set(FFTW_LIB              "/usr/lib64/libfftw3.so")
set(NETCDF_INCLUDE_DIR    "/home/nyue/Systems/netcdf/4.6.1/include")
set(NETCDFCXX_INCLUDE_DIR "/home/nyue/Systems/netcdf-cxx/4.3.0/include")
set(NETCDF_LIB_C          "/home/nyue/Systems/netcdf/4.6.1/lib64/libnetcdf.so")
set(NETCDF_LIB_CPP        "/home/nyue/Systems/netcdf-cxx/4.3.0/lib64/libnetcdf-cxx4.so")
set(HDF5_LIB_1            "/usr/lib64/libhdf5.so")
set(HDF5_LIB_2            "/usr/lib64/libhdf5_hl.so")
set(SZIP_LIB              "/usr/lib64/libsz.so")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR} ${NETCDFCXX_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=__restrict__)
