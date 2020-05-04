# LISA GPU nodes
# Tested with:
# module load CUDA
# module load Boost/1.65.1-foss-2017b
# module load netCDF/4.5.0-foss-2017b 
# module load Python/3.6.3-foss-2017b
# module load CMake/3.10.2-GCCcore-6.4.0
# module unload ScaLAPACK/2.0.2-gompi-2017b-OpenBLAS-0.2.20
# module unload gompi/2017b

set(ENV{CC}  gcc) # C compiler for serial build
set(ENV{CXX} g++) # C++ compiler for serial build

# set(USER_CXX_FLAGS " -std=c++14 -D_GLIBCXX_USE_CXX11_ABI=0")
set(USER_CXX_FLAGS "-std=c++14 -fopenmp")
set(USER_CXX_FLAGS_RELEASE "-Ofast -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_LIB     "fftw3")
set(FFTWF_LIB    "fftw3f")
set(NETCDF_LIB_C "netcdf")
set(HDF5_LIB     "hdf5")
set(SZIP_LIB     "sz")

set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB} ${SZIP_LIB} m z curl)

if(USECUDA)
   set(CUDA_PROPAGATE_HOST_FLAGS OFF)
   # set(LIBS ${LIBS} -rdynamic /hpc/sw/cuda/8.0.44/lib64/libcufft.so)
   set(LIBS ${LIBS} -rdynamic libcufft.so)
   set(USER_CUDA_NVCC_FLAGS "-arch=sm_70")
  list(APPEND CUDA_NVCC_FLAGS " -std=c++14")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
