# SuperMUC
# set(ENV{CXX} icc) # compiler for serial build
set(ENV{CXX} mpCC) # compiler for parallel build

set(USER_CXX_FLAGS "-restrict -DMPICH_IGNORE_CXX_SEEK")
set(USER_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -xHOST -fno-alias -restrict -vec-report1 -no-prec-div")
set(USER_CXX_FLAGS_DEBUG "-traceback -check=conversions,stack,uninit -check-pointers=rw -check-pointers-dangling=all-check-pointers-undimensioned -fp-stack-check -fp-trap=common -fp-trap-all=common")

set(FFTW_INCLUDE_DIR   "/lrz/sys/libraries/fftw/3.3.2/sse/include")
set(FFTW_LIB           "/lrz/sys/libraries/fftw/3.3.2/sse/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/lrz/sys/libraries/netcdf/4.2.1.1/include")
set(NETCDF_LIB_C       "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5.a")
set(HDF5_LIB_2         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5_hl.a")
set(SZIP_LIB           "/lrz/sys/libraries/hdf5/szip_2.1_u1/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=restrict)
