# MPI-PC, Intel compilers
if(USEMPI) 
  set(ENV{CC}  mpicc) # C compiler for parallel build
  set(ENV{CXX} mpiCC) # C++ compiler for parallel build
else()
  set(ENV{CC}  icc ) # C compiler for serial build
  set(ENV{CXX} icpc) # C++ compiler for serial build
endif()

set(USER_CXX_FLAGS "")
set(USER_CXX_FLAGS "-restrict")
set(USER_CXX_FLAGS_RELEASE "-xHOST -O3")
set(USER_CXX_FLAGS_DEBUG "-traceback -check=conversions,stack,uninit -check-pointers=rw -check-pointers-dangling=all-check-pointers-undimensioned -fp-stack-check -fp-trap=common -fp-trap-all=common")

set(FFTW_INCLUDE_DIR   "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/include")
set(FFTW_LIB           "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/sw/squeeze-x64/netcdf-latest-static-intel13/include")
set(NETCDF_LIB_C       "/sw/squeeze-x64/netcdf-latest-static-intel13/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/sw/squeeze-x64/netcdf-latest-static-intel13/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5.a")
set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5_hl.a")
set(SZIP_LIB           "/sw/squeeze-x64/szip-latest-static/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=restrict)
