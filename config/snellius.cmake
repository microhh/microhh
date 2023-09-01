# Snellius @ SURFSara.

# This is an example script that load the correct modules from 2022 stack.
##!/bin/sh
# module purge
# module load 2022
# module load CMake/3.23.1-GCCcore-11.3.0
# 
# # GCC
# module load foss/2022a
# module load netCDF/4.9.0-gompi-2022a
# module load CUDA/11.8.0
# module load Clang/13.0.1-GCCcore-11.3.0
# 
# # Python et al.
# module load ncview/2.1.8-gompi-2022a
# module load Python/3.10.4-GCCcore-11.3.0
# module load IPython/8.5.0-GCCcore-11.3.0
# module load NCO/5.1.0-foss-2022a
# module load Tk/8.6.12-GCCcore-11.3.0
# End example script


# Select correct compilers for GCC + parallel/serial:
if(USEMPI)
    set(ENV{CXX} mpicxx)
    set(ENV{FC} mpif90)
else()
    set(ENV{CXX} g++)
    set(ENV{FC} gfortran)
endif()

# Set compiler flags / options:
if(USECUDA)
    set(USER_CXX_FLAGS "-fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-O3 -march=icelake-server -mtune=icelake-server")
    set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    set(USER_CXX_FLAGS "")
    set(USER_CXX_FLAGS_RELEASE "-O3 -march=znver2 -mtune=znver2 -mfma -mavx2 -m3dnow -fomit-frame-pointer")
    set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

    set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
    set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=znver2 -mtune=znver2 -mfma -mavx2 -m3dnow -fomit-frame-pointer")

    add_definitions(-DRESTRICTKEYWORD=__restrict__)
endif()


set(NETCDF_LIB_C "netcdf")
set(FFTW_LIB "fftw3")
set(FFTWF_LIB "fftw3f")
set(IRC_LIB "irc")
set(IRC_LIB "")
set(HDF5_LIB "hdf5")
set(SZIP_LIB "sz")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB} ${SZIP_LIB} ${IRC_LIB} m z curl)


if(USECUDA)
    set(CMAKE_CUDA_ARCHITECTURES 80)
    set(USER_CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
    set(USER_CUDA_NVCC_FLAGS_RELEASE "-Xptxas -O3 -DNDEBUG")
    set(USER_CUDA_NVCC_FLAGS_DEBUG "-Xptxas -O0 -g -DCUDACHECKS")
    add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_CUDA)
endif()

# Disable MPI-IO for cross-sections on GPFS file systems.
add_definitions(-DDISABLE_2D_MPIIO=1)
add_definitions(-DRTE_USE_CBOOL)
