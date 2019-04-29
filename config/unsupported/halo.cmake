# HALO, DKRZ Visualization cluster
if(USEMPI) 
  set(ENV{CC}  /sw/centos58-x64/mpi/openmpi-1.8.1-cuda-shared-gcc47/bin/mpicc) # C compiler for serial build
  set(ENV{CXX} /sw/centos58-x64/mpi/openmpi-1.8.1-cuda-shared-gcc47/bin/mpicxx) # C++ compiler for serial build
else()
  set(ENV{CC}  gcc) # C compiler for serial build
  set(ENV{CXX} g++) # C++ compiler for serial build
endif()

set(USER_C_FLAGS "")
set(USER_C_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_C_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(USER_CXX_FLAGS "")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/home/zmaw/m300207/include")
set(FFTW_LIB           "/home/zmaw/m300207/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/sw/centos58-x64/netcdf-4.1.3-static-gcc47/include")
set(NETCDF_LIB_C       "/sw/centos58-x64/netcdf-4.1.3-static-gcc47/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/sw/centos58-x64/netcdf-4.1.3-static-gcc47/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/sw/centos58-x64/hdf5-latest-static/lib/libhdf5.a")
set(HDF5_LIB_2         "/sw/centos58-x64/hdf5-latest-static/lib/libhdf5_hl.a")
set(SZIP_LIB           "/sw/centos58-x64/szip-latest-static/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(LIBS ${LIBS} -rdynamic /usr/local/cuda/lib64/libcufft.so)
  set(USER_CUDA_NVCC_FLAGS "-arch=sm_20")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
