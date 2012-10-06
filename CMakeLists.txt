cmake_minimum_required (VERSION 2.6)
project (microhhc)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# MPI PC
set(FFTW_INCLUDE_DIR   "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/include")
set(FFTW_LIB           "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/sw/squeeze-x64/netcdf-latest-static-gcc46/include")
set(NETCDF_LIB_C       "/sw/squeeze-x64/netcdf-latest-static-gcc46/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/sw/squeeze-x64/netcdf-latest-static-gcc46/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5.a")
set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5_hl.a")
set(SZIP_LIB           "/sw/squeeze-x64/szip-latest-static/lib/libsz.a")

# Macbook
set(FFTW_INCLUDE_DIR   "/opt/local/include")
set(FFTW_LIB           "/opt/local/lib/libfftw3.dylib")
set(NETCDF_INCLUDE_DIR "/opt/local/include")
set(NETCDF_LIB_C       "/opt/local/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/opt/local/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/opt/local/lib/libhdf5.dylib")
set(HDF5_LIB_2         "/opt/local/lib/libhdf5_hl.dylib")
set(SZIP_LIB           "/opt/local/lib/libsz.a")

#SuperMUC
#set(FFTW_INCLUDE_DIR   "/lrz/sys/libraries/fftw/3.3.2/sse/include")
#set(FFTW_LIB           "/lrz/sys/libraries/fftw/3.3.2/sse/lib/libfftw3.a")
#set(NETCDF_INCLUDE_DIR "/lrz/sys/libraries/netcdf/4.2.1.1/include")
#set(NETCDF_LIB_C       "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf.a")
#set(NETCDF_LIB_CPP     "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf_c++.a")
#set(HDF5_LIB_1         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5.a")
#set(HDF5_LIB_2         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5_hl.a")
#set(SZIP_LIB           "/lrz/sys/libraries/hdf5/szip_2.1_u1/lib/libsz.a")

#set(FFTW_INCLUDE_DIR "/bgsys/local/fftw3/3.3/fftw_450/include")
#set(FFTW_LIB         "/bgsys/local/fftw3/3.3/fftw_450/lib/libfftw3.a")

#set(FFTW_INCLUDE_DIR "/bgsys/local/fftw3/3.3.2/fftw_g/include")
#set(FFTW_LIB         "/bgsys/local/fftw3/3.3.2/fftw_g/lib/libfftw3.a")

set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB})

add_subdirectory(src)
add_subdirectory(main)
add_subdirectory(test)

