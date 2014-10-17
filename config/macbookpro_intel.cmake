# MPI-PC
#set(ENV{CXX} icpc) # compiler for serial build
set(ENV{CXX} mpicxx) # compiler for parallel build
set(GNU_SED "gsed")

set(CMAKE_BUILD_TYPE "Release")
set(USER_CXX_FLAGS "-restrict") 
set(USER_CXX_FLAGS_RELEASE "-xHOST -O3")
#set(USER_CXX_FLAGS_DEBUG "-traceback -check=conversions,stack,uninit -check-pointers=rw -check-pointers-dangling=all-check-pointers-undimensioned -fp-stack-check -fp-trap=common -fp-trap-all=common")
set(USER_CXX_FLAGS_DEBUG "-debug -g -check=conversions,stack,uninit -fp-stack-check -fp-trap=common -fp-trap-all=common ") 
#"-traceback -check=conversions,stack,uninit -fp-stack-check -fp-trap=common -fp-trap-all=common")

set(FFTW_INCLUDE_DIR   "/usr/local/include")
set(FFTW_LIB           "/usr/local/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/usr/local/netcdf/include/")
set(NETCDF_LIB_C       "/usr/local/netcdf/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/usr/local/netcdf/lib/libnetcdf_c++.a")
#set(HDF5_LIB_1         "/usr/local/hdf/lib/libhdf5.a")
#set(HDF5_LIB_2         "/usr/local/hdf/lib/libhdf5_hl.a")
#set(SZIP_LIB           "/usr/local/szip/lib/libsz.a")
#set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)

set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} m z curl)

add_definitions(-DRESTRICTKEYWORD=restrict)
