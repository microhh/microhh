# Blizzard
# set(ENV{CXX} xlc++) # compiler for serial build
set(ENV{CXX} mpCC) # compiler for parallel build

set(USER_CXX_FLAGS_RELEASE "-qarch=pwr6 -qtune=pwr6 -O3 -qhot=simd -qenablevmx")

set(FFTW_INCLUDE_DIR   "/pf/zmaw/m300041/local/include")
set(FFTW_LIB           "/pf/zmaw/m300041/local/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/sw/aix61/netcdf-4.1.3/include")
set(NETCDF_LIB_C       "/sw/aix61/netcdf-4.1.3/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/sw/aix61/netcdf-4.1.3/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/sw/aix61/hdf5-1.8.8/lib/libhdf5.a")
set(HDF5_LIB_2         "/sw/aix61/hdf5-1.8.8/lib/libhdf5_hl.a")
set(SZIP_LIB           "/sw/aix53/szip-2.1/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=__restrict__)
