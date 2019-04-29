# Mistral GPU nodes
if(USEMPI) 
  set(ENV{CC}  mpiicc ) # C compiler for parallel build
  set(ENV{CXX} mpiicpc) # C++ compiler for parallel build
else()
  set(ENV{CC}  gcc) # C compiler for serial build
  set(ENV{CXX} g++) # C++ compiler for serial build
endif()

set(USER_CXX_FLAGS "-restrict")
set(USER_CXX_FLAGS_RELEASE "-Ofast -xCORE-AVX2")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

# Note(s) BvS: * the legacy NetCDF C++ library is not (yet) installed on Mistral
#              * FFTW3 with NVCC suffers from this problem: https://github.com/FFTW/fftw3/issues/18
# To solve both issues, config file uses libraries from my home directory
set(FFTW_INCLUDE_DIR   "/home/zmaw/m300241/tools/fftw3/include")
set(FFTW_LIB           "/home/zmaw/m300241/tools/fftw3/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR1 "/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-static-gcc48/include")
set(NETCDF_INCLUDE_DIR2 "/home/zmaw/m300241/tools/netcdfcpp/include")
set(NETCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR1} ${NETCDF_INCLUDE_DIR2})
set(NETCDF_LIB_C       "/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-static-gcc48/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/home/zmaw/m300241/tools/netcdfcpp/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/sw/rhel6-x64/hdf5/hdf5-1.8.14-gcc48/lib/libhdf5.a")
set(HDF5_LIB_2         "/sw/rhel6-x64/hdf5/hdf5-1.8.14-gcc48/lib/libhdf5_hl.a")
set(SZIP_LIB           "/sw/rhel6-x64/sys/szip-2.1-gcc48/lib/libsz.a")
set(CURL_LIB           "/usr/lib64/libcurl.so.4")
set(LIBS ${CURL_LIB} ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z dl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(LIBS ${LIBS} -rdynamic /opt/cuda/6.5/lib64/libcufft.so)
  set(USER_CUDA_NVCC_FLAGS "-arch=sm_37")
endif()

add_definitions(-DRESTRICTKEYWORD=restrict)
