# OSC-Owens cluster
if(USECUDA)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(CUDALIBS "-rdynamic /usr/local/cuda/9.1.85/lib64/libcufft.so")
  set(USER_CUDA_NVCC_FLAGS "--expt-relaxed-constexpr -arch=sm_60")
  add_definitions(-DRESTRICTKEYWORD=__restrict__)
  list(APPEND CUDA_NVCC_FLAGS "-std=c++14")


set(USER_C_FLAGS " -std=c++14")
set(USER_C_FLAGS_RELEASE "-Ofast -march=native")
set(USER_C_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(USER_CXX_FLAGS " -std=c++14")
set(USER_CXX_FLAGS_RELEASE "-Ofast -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
else()
  add_definitions(-DRESTRICTKEYWORD=restrict)



set(USER_C_FLAGS "--expt-relaxed-constexpr -std=c++14")
set(USER_C_FLAGS_RELEASE "-Ofast -march=native")
set(USER_C_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")


set(USER_CXX_FLAGS "--expt-relaxed-constexpr  -std=c++14")
set(USER_CXX_FLAGS_RELEASE "-Ofast -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
endif()

set(USER_FC_FLAGS "-debug -traceback -r8 -ftz -extend_source ")
set(USER_FC_FLAGS_RELEASE " -O3 -no-prec-div  -xHost")
set(USER_FC_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created ")

if(USEMPI)
  set(ENV{CC}  mpicc) # C compiler for serial build
  set(ENV{CXX} mpicc) # C++ compiler for serial build
  set(ENV{FC}  mpif90) # Fortran compiler for serial build

  set(FFTW_INCLUDE_DIR   "/usr/local/fftw3/intel/17.0/mvapich2/2.3rc2/3.3.5/include/")
  set(FFTW_LIB           "/usr/local/fftw3/intel/17.0/mvapich2/2.3rc2/3.3.5/lib/libfftw3.a")
  set(FFTWF_LIB          "/usr/local/fftw3/intel/17.0/mvapich2/2.3rc2/3.3.5/lib/libfftw3f.a")
  set(NETCDF_INCLUDE_DIR "/usr/local/netcdf/intel/17.0/mvapich2/2.3rc2/4.3.3.1/include/")
  set(NETCDF_LIB_C       "/usr/local/netcdf/intel/17.0/mvapich2/2.3rc2/4.3.3.1/lib/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/usr/local/netcdf/intel/17.0/mvapich2/2.3rc2/4.3.3.1/lib/libnetcdf_c++4.a")
  set(HDF5_LIB_1         "/usr/local/hdf5/intel/17.0/mvapich2/2.3rc2/1.8.17/lib/libhdf5.a")
  set(HDF5_LIB_2         "/usr/local/hdf5/intel/17.0/mvapich2/2.3rc2/1.8.17/lib/libhdf5_hl.a")
  set(SZIP_LIB           "")
else()
  set(ENV{CC} gcc) # compiler for serial build
  set(ENV{CXX} g++) # compiler for parallel build
  set(ENV{FC}  gfortran)
  set(FFTW_INCLUDE_DIR   "/usr/local/fftw3/gnu/6.3/mvapich2/2.3rc2/3.3.5/include/")
  set(FFTW_LIB           "/usr/local/fftw3/gnu/6.3/mvapich2/2.3rc2/3.3.5/lib/libfftw3.a")
  set(FFTWF_LIB          "/usr/local/fftw3/gnu/6.3/mvapich2/2.3rc2/3.3.5/lib/libfftw3f.a")
   set(NETCDF_INCLUDE_DIR "/usr/local/netcdf/gnu/6.3/4.3.3.1-serial/include")
  set(NETCDF_LIB_C       "/usr/local/netcdf/gnu/6.3/4.3.3.1-serial/lib/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/usr/local/netcdf/gnu/6.3/4.3.3.1-serial/lib/libnetcdf_c++4.a")
  set(HDF5_LIB_1         "/usr/local/hdf5/gnu/6.3/1.8.17-serial/lib/libhdf5.a")
  set(HDF5_LIB_2         "/usr/local/hdf5/gnu/6.3/1.8.17-serial/lib/libhdf5_hl.a")
  set(SZIP_LIB           "")
endif()

set(LIBS ${CUDALIBS} ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR} ${NETCDF_INCLUDE_CXX_DIR})
