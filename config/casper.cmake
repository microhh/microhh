# NCAR Casper, GPU-enabled, gnu compilers
#
# execdav --constraint=gp100
# module load gnu/8.3.0 
# module load netcdf/4.7.3
# module load fftw/3.3.8
# module load cuda/10.1
# 
# note: Casper does not provide fftw with single precision support
#       
set(ENV{CC}  gcc ) # C compiler for serial build
set(ENV{CXX} g++ ) # C++ compiler for serial build

set(USER_CXX_FLAGS "-std=c++14")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-traceback ")
set(BOOST_INCLUDE_DIR  "/glade/u/apps/ch/opt/boost/1_67_0/include")
set(FFTW_INCLUDE_DIR   "/glade/u/apps/dav/opt/fftw/3.3.8/gnu/8.3.0/include")
set(FFTW_LIB           "/glade/u/apps/dav/opt/fftw/3.3.8/gnu/8.3.0/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/glade/u/apps/dav/opt/netcdf/4.7.3/gnu/8.3.0/include")
set(NETCDF_LIB_C       "/glade/u/apps/dav/opt/netcdf/4.7.3/gnu/8.3.0/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/glade/u/apps/dav/opt/netcdf/4.7.3/gnu/8.3.0/lib/libnetcdf_c++4.a")
set(HDF5_LIB_1         "/glade/u/apps/dav/opt/netcdf/4.7.3/gnu/8.3.0/lib/libhdf5.a")
set(HDF5_LIB_2         "/glade/u/apps/dav/opt/netcdf/4.7.3/gnu/8.3.0/lib/libhdf5_hl.a")
set(SZIP_LIB           "/glade/u/apps/dav/opt/netcdf/4.7.3/gnu/8.3.0/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${BOOST_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(CUFFT_LIB "/glade/u/apps/dav/opt/cuda/10.1/lib64/libcufft.so")
  set(LIBS ${LIBS} ${CUFFT_LIB} -rdynamic )
  set(USER_CUDA_NVCC_FLAGS "-arch=sm_70")
  list(APPEND CUDA_NVCC_FLAGS "-std=c++14")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
