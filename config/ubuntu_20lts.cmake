# CMake file for Ubuntu 20 and 22 LTS based on generic.cmake.

if(USEMPI)
    set(ENV{CC} mpicc )
    set(ENV{CXX} mpicxx)
    set(ENV{FC} mpif90)
else()
    set(ENV{CC} gcc )
    set(ENV{CXX} g++)
    set(ENV{FC} gfortran)
endif()

# Set compiler flags / options:
if(USECUDA)
    set(USER_CXX_FLAGS "-std=c++17 -fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-O3")
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    set(USER_CXX_FLAGS "-std=c++17")
    set(USER_CXX_FLAGS_RELEASE "-O3 -march=native")

    set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
    set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")

    add_definitions(-DRESTRICTKEYWORD=__restrict__)
endif()

set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(NETCDF_LIB_C "netcdf")
set(FFTW_LIB "fftw3")
set(FFTWF_LIB "fftw3f")
set(HDF5_LIB "hdf5_serial")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB}) #It may be necessary to add m z curl sz if necessary

if(USECUDA)
    set(CMAKE_CUDA_ARCHITECTURES 86)
    set(USER_CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
    set(USER_CUDA_NVCC_FLAGS_RELEASE "-DNDEBUG")
    set(USER_CUDA_NVCC_FLAGS_DEBUG "-O0 -g -DCUDACHECKS")
    add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_CUDA)
endif()

# Disable MPI-IO for cross-sections on GPFS file systems. This may or may not be necessary, depending on the system
add_definitions(-DDISABLE_2D_MPIIO=1)
add_definitions(-DRTE_USE_CBOOL)
