
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
    set(USER_CXX_FLAGS "-std=c++17 -fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-O3 -march=icelake-server -mtune=icelake-server")
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    if(USEINTEL)
        set(USER_CXX_FLAGS "-std=c++17 -restrict")
        set(USER_CXX_FLAGS_RELEASE "-O3 -march=core-avx2")
        add_definitions(-DRESTRICTKEYWORD=restrict)
    else()
        set(USER_CXX_FLAGS "-std=c++17")
        set(USER_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
        set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
        set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
        set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")
        set(USER_FC_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

        add_definitions(-DRESTRICTKEYWORD=__restrict__)
    endif()
endif()

# Set compiler flags / options:
if(USECUDA)
    set(USER_CXX_FLAGS "-std=c++17 -fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-O3 -march=icelake-server -mtune=icelake-server")
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    if(USEINTEL)
        set(USER_CXX_FLAGS "-std=c++17 -restrict")
        set(USER_CXX_FLAGS_RELEASE "-O3 -march=core-avx2")
        add_definitions(-DRESTRICTKEYWORD=restrict)
    else()
        set(USER_CXX_FLAGS "-std=c++17")
        set(USER_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
        set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
        set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
        set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")
        set(USER_FC_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

        add_definitions(-DRESTRICTKEYWORD=__restrict__)
    endif()
endif()


set(BOOST_INCLUDE_DIR  "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/boost-1.77.0-gvqdfio5uj6a625h64y3jppv3jjfzx7m/include")

set(NETCDF_LIB_C "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/netcdf-c-4.8.1-kz7m3osaphp3uut6i2tg5a5mdqf7q64m/lib/libnetcdf.so")
set(NETCDF_INCLUDE_DIR "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/netcdf-c-4.8.1-kz7m3osaphp3uut6i2tg5a5mdqf7q64m/include")


set(FFTW_LIB "/apps/noarch/software/fftw3/lib/libfftw3.so")
#set(FFTWF_LIB "fftw3f")

set(FFTW_INCLUDE_DIR "/apps/noarch/software/fftw3/include")


set(IRC_LIB "irc")
set(IRC_LIB "")

set(HDF5_LIB_1 "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/hdf5-1.10.7-wscpmjfq75bppp3geu4xtecw3buxhnke/lib/libhdf5.so")
set(HDF5_LIB_2 "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/hdf5-1.10.7-wscpmjfq75bppp3geu4xtecw3buxhnke/lib/libhdf5_hl.so")

set(SZIP_LIB "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/sz-2.1.12-lugn4xqq2qgwmtl3g3irziro4zazyubk/lib64/libSZ.so")
set(SZIP_INCLUDE_DIR "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/sz-2.1.12-lugn4xqq2qgwmtl3g3irziro4zazyubk/include")

set(CURL_LIB "/apps/arch/2022r2/software/linux-rhel8-skylake_avx512/gcc-8.5.0/curl-7.79.0-kirivlbuvyhlieiv6qj67yyl4lm36cbo/lib/libcurl.so")

set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} ${IRC_LIB} ${CURL_LIB} m z curl)
set(INCLUDE_DIRS ${BOOST_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR} ${SZIP_INCLUDE_DIR})
    
# Disable MPI-IO for cross-sections on GPFS file systems.
# add_definitions(-DDISABLE_2D_MPIIO=1)

if(USECUDA)
    set(CMAKE_CUDA_ARCHITECTURES 80)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    set(LIBS ${LIBS} -rdynamic $ENV{EBROOTCUDA}/lib64/libcufft.so)
    set(USER_CUDA_NVCC_FLAGS "-std=c++17 --expt-relaxed-constexpr")
    set(USER_CUDA_NVCC_RELEASE_FLAGS "-O3 --use_fast_math")
    add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_CUDA)
endif()

add_definitions(-DRTE_RRTMGP_USE_CBOOL)
