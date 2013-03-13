# ARCH Linux
set(CMAKE_C_COMPILER   "/usr/bin/cc")
set(CMAKE_CXX_COMPILER "/usr/bin/c++")
set (CXX_COMPILER_WRAPPER mpicxx)
set (C_COMPILER_WRAPPER mpicc)

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(FFTW_INCLUDE_DIR   "/usr/include")
set(FFTW_LIB           "/usr/lib64/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_C       "/usr/lib64/libnetcdf.a")
set(NETCDF_LIB_CPP     "/usr/lib64/libnetcdf_c++.a")
set(HDF5_LIB_1         "/usr/lib64/libhdf5.a")
set(HDF5_LIB_2         "/usr/lib64/libhdf5_hl.a")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
