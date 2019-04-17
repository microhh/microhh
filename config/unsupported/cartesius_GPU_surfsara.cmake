# Cartesius GPU nodes
# Tested with:
#module load compilerwrappers
#module load c/intel/15.0.0
#module load fortran/intel/15.0.0
#module load mpi/impi/5.0.3.048
#module load bull
#module load surfsara
#module load gcc/5.2.0
#module load cuda/8.0.44 
#module load hdf5/serial/intel/1.8.10-patch1
#module load netcdf-cxx/serial/intel/4.3.0
#module load netcdf/serial/intel/4.3.3.1 
#module load szip/gnu/2.1
#module load fftw3/intel/3.3.3

set(ENV{CC}  gcc) # C compiler for serial build
set(ENV{CXX} g++) # C++ compiler for serial build

set(USER_CXX_FLAGS " -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0")
set(USER_CXX_FLAGS_RELEASE "-Ofast -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

#set(FFTW_INCLUDE_DIR "/home/bstratum/tools/fftw3_linked/include")
#set(FFTW_LIB         "/home/bstratum/tools/fftw3_linked/lib/libfftw3.a")
#set(HDF5_LIB_1 "/home/bstratum/tools/hdf5-1.8.17-gcc480/lib/libhdf5.a")
#set(HDF5_LIB_2 "/home/bstratum/tools/hdf5-1.8.17-gcc480/lib/libhdf5_hl.a")

set(FFTW_LIB       "fftw3")
set(FFTWF_LIB      "fftw3f")
set(NETCDF_LIB_C   "netcdf")
set(NETCDF_LIB_CPP "netcdf_c++4")
set(IRC_LIB        "irc")
set(HDF5_LIB       "hdf5")
set(SZIP_LIB       "sz")

set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB} ${SZIP_LIB} ${IRC_LIB} m z curl)
#set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR} ${NETCDF_INCLUDE_CXX_DIR})

if(USECUDA)
   set(CUDA_PROPAGATE_HOST_FLAGS OFF)
   set(LIBS ${LIBS} -rdynamic /hpc/sw/cuda/8.0.44/lib64/libcufft.so)
   set(USER_CUDA_NVCC_FLAGS "-arch=sm_35")
  list(APPEND CUDA_NVCC_FLAGS "-std=c++11")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
