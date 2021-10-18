if(USEMPI)
    set(ENV{CC}  mpiicc ) # C compiler for parallel build
    set(ENV{CXX} mpiicpc) # C++ compiler for parallel build
    set(USEICC TRUE)
else()
    set(ENV{CC}  icc ) # C compiler for parallel build
    set(ENV{CXX} icpc) # C++ compiler for serial build
    set(USEICC TRUE)
endif()

if(USECUDA)
    set(ENV{CC}  gcc) # C compiler for serial build
    set(ENV{CXX} g++) # C++ compiler for serial build
    set(USEICC FALSE)

    set(USER_CXX_FLAGS "-std=c++14 -fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-Ofast -march=ivybridge") # -march optimized for the CPU present in Cartesius GPU nodes
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    set(USER_CXX_FLAGS "-std=c++14 -restrict")
    set(USER_CXX_FLAGS_RELEASE "-Ofast")
    add_definitions(-DRESTRICTKEYWORD=restrict)
endif()

set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_LIB "fftw3")
set(FFTWF_LIB "fftw3f")
set(NETCDF_LIB_C "netcdf")
set(IRC_LIB "irc")
set(IRC_LIB "")
set(HDF5_LIB "hdf5")
set(SZIP_LIB "sz")

set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB} ${SZIP_LIB} ${IRC_LIB} m z curl)

if(USECUDA)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    set(LIBS ${LIBS} -rdynamic $ENV{EBROOTCUDA}/lib64/libcufft.so)
    set(USER_CUDA_NVCC_FLAGS "-arch=sm_80")
    list(APPEND CUDA_NVCC_FLAGS "-std=c++14")
    list(APPEND CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
endif()
