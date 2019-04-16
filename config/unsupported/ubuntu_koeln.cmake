# ARCH Linux
set(CMAKE_C_COMPILER   "/usr/bin/cc")
set(CMAKE_CXX_COMPILER "/usr/bin/c++")
set(CXX_COMPILER_WRAPPER mpicxx)
set(C_COMPILER_WRAPPER mpicc)

set(USER_CXX_FLAGS "-fpermissive ")
set(USER_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/usr/local/include")
set(NETCDF_LIB_C       "/usr/local/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/usr/local/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/usr/lib/libhdf5.a")
set(HDF5_LIB_2         "/usr/lib/libhdf5_hl.a")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=__restrict__)
