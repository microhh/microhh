# ARCH Linux / Manjaro
if(USEMPI)
  set(ENV{CC}  mpicc)    # C compiler for parallel build
  set(ENV{CXX} mpicxx)   # C++ compiler for parallel build
  set(ENV{FC}  mpif90)   # Fortran compiler for parallel build
else()
  if(USECUDA)
    set(ENV{CC}  gcc-11)      # C compiler for serial build
    set(ENV{CXX} g++-11)      # C++ compiler for serial build
    set(ENV{FC}  gfortran-11) # Fortran compiler for serial build
  else()
    set(ENV{CC}  gcc)        # C compiler for serial build
    set(ENV{CXX} g++)        # C++ compiler for serial build
    set(ENV{FC}  gfortran)   # Fortran compiler for serial build
  endif()
endif()

if(USECUDA)
  set(USER_CXX_FLAGS "-std=c++17 -fopenmp")
else()
  set(USER_CXX_FLAGS "-std=c++17")
endif()

set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
set(USER_FC_FLAGS_RELEASE "-DNDEBUG -Ofast -march=native")
set(USER_FC_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib/libfftw3.so")
set(FFTWF_LIB          "/usr/lib/libfftw3f.so")
set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_C       "/usr/lib/libnetcdf.so")
set(HDF5_LIB_1         "/usr/lib/libhdf5.so")
set(HDF5_LIB_2         "/usr/lib/libhdf5_hl.so")
set(SZIP_LIB           "/usr/lib/libsz.so")

set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl dl rt)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

if(USECUDA)
    set(CMAKE_CUDA_ARCHITECTURES 70)
    set(USER_CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
    set(USER_CUDA_NVCC_FLAGS_RELEASE "-DNDEBUG")
    set(USER_CUDA_NVCC_FLAGS_DEBUG "-O0 -g -DCUDACHECKS")
    add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_CUDA)
endif()

add_definitions(-DRTE_USE_CBOOL)
add_definitions(-DRESTRICTKEYWORD=__restrict__)
