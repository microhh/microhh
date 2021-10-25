# Switch between Intel and GCC:
set(USEINTEL FALSE)

# GPU builds are always with GCC:
if(USECUDA)
    set(USEINTEL FALSE)
endif()

# Select correct compilers for Intel/GCC + parallel/serial:
if(USEMPI)
    if(USEINTEL)
        set(ENV{CC} mpiicc )
        set(ENV{CXX} mpiicpc)
	set(ENV{FC} mpiifort)
    else()
        set(ENV{CC} mpicc )
        set(ENV{CXX} mpicxx)
        set(ENV{FC} mpif90)
    endif()
else()
    if(USEINTEL)
        set(ENV{CC} icc )
        set(ENV{CXX} icpc)
        set(ENV{FC} ifort)
    else()
        set(ENV{CC} gcc )
        set(ENV{CXX} g++)
        set(ENV{FC} gfortran)
    endif()
endif()

# Set compiler flags / options:
if(USECUDA)
    set(USER_CXX_FLAGS "-std=c++14 -fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-Ofast -march=icelake-server -mtune=icelake-server")
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    if(USEINTEL)
        set(USER_CXX_FLAGS "-std=c++14 -restrict")
        set(USER_CXX_FLAGS_RELEASE "-Ofast -march=core-avx2")
        add_definitions(-DRESTRICTKEYWORD=restrict)
    else()
        set(USER_CXX_FLAGS "-std=c++14")
        set(USER_CXX_FLAGS_RELEASE "-Ofast -march=znver2 -mtune=znver2 -mfma -mavx2 -m3dnow -fomit-frame-pointer")

        set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
        set(USER_FC_FLAGS_RELEASE "-DNDEBUG -Ofast -march=znver2 -mtune=znver2 -mfma -mavx2 -m3dnow -fomit-frame-pointer")

        add_definitions(-DRESTRICTKEYWORD=__restrict__)
    endif()
endif()

set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

if(USEINTEL)
    set(NETCDF_INCLUDE_DIR "/home/bstratum/software/netcdf-c-4.8.1-intel-2021a/include")
    set(NETCDF_LIB_C "/home/bstratum/software/netcdf-c-4.8.1-intel-2021a/lib/libnetcdf.so")
    set(INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
else()
    set(NETCDF_LIB_C "netcdf")
endif()

set(FFTW_LIB "fftw3")
set(FFTWF_LIB "fftw3f")
set(IRC_LIB "irc")
set(IRC_LIB "")
set(HDF5_LIB "hdf5")
set(SZIP_LIB "sz")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB} ${SZIP_LIB} ${IRC_LIB} m z curl)

# Disable MPI-IO for cross-sections on GPFS file systems.
add_definitions(-DDISABLE_2D_MPIIO=1)

if(USECUDA)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    set(LIBS ${LIBS} -rdynamic $ENV{EBROOTCUDA}/lib64/libcufft.so)
    set(USER_CUDA_NVCC_FLAGS "-arch=sm_80 --use_fast_math")
    list(APPEND CUDA_NVCC_FLAGS "-std=c++14 --expt-relaxed-constexpr")
endif()
