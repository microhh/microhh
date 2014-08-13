# MPI-PC
# set(ENV{CXX} g++) # compiler for serial build
set(ENV{CXX} /sw/centos58-x64/mpi/openmpi-1.8.1-cuda-shared-gcc47/bin/mpicxx) # compiler for parallel build

set(USER_CXX_FLAGS "")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/sw/centos58-x64/numerics/fftw-3.3-openmp-gcc46/include")
set(FFTW_LIB           "/sw/centos58-x64/numerics/fftw-3.3-openmp-gcc46/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/sw/centos58-x64/netcdf-4.1.3-static-gcc47/include")
set(NETCDF_LIB_C       "/sw/centos58-x64/netcdf-4.1.3-static-gcc47/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/sw/centos58-x64/netcdf-4.1.3-static-gcc47/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/sw/centos58-x64/hdf5-latest-static/lib/libhdf5.a")
set(HDF5_LIB_2         "/sw/centos58-x64/hdf5-latest-static/lib/libhdf5_hl.a")
set(SZIP_LIB           "/sw/centos58-x64/szip-latest-static/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)

if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(CUDA_NVCC_FLAGS "-arch=sm_21")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
