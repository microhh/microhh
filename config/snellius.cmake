#
# Snellius @ SURFSara.
# For instructions, best practices, and an overview of the required software modules, see:
# https://microhh.readthedocs.io/en/latest/computing_systems/snellius.html
#

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
