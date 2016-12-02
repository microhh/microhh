# WUR HPC cluster
if(USEMPI) 
  set(ENV{CXX} mpiicpc) # compiler for parallel build
else()
  set(ENV{CXX} icpc) # C++ compiler for serial build
endif()

set(USER_CXX_FLAGS "-restrict -std=c++11")
set(USER_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -xHOST -fno-alias")
set(USER_CXX_FLAGS_DEBUG "-traceback -check=conversions,stack,uninit -check-pointers=rw -check-pointers-dangling=all-check-pointers-undimensioned -fp-stack-check -fp-trap=common -fp-trap-all=common")

set(FFTW_INCLUDE_DIR    "/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include")
set(FFTW_LIB            "/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR  "/cm/shared/apps/netcdf/intel/64/4.4.1/include")
set(NETCDF_LIB_C        "/cm/shared/apps/netcdf/intel/64/4.4.1/lib/libnetcdf.a")
set(NETCDF_LIB_CPP      "/cm/shared/apps/netcdf/intel/64/4.4.1/lib/libnetcdf_c++4.a")
set(HDF5_LIB_1          "/cm/shared/apps/hdf5/intel/64/1.8.17/lib/libhdf5.a")
set(HDF5_LIB_2          "/cm/shared/apps/hdf5/intel/64/1.8.17/lib/libhdf5_hl.a")
set(SZIP_LIB            "/cm/shared/apps/szip/intel/64/2.1/lib/libsz.a")

set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=restrict)
