# ECMWF CCA/CCB
# prgenvswitchto intel
# module load cmake
# module load netcdf4
# module load fftw
# module load szip
# module load cray-hdf5
# module load gcc
# module load boost

if(USEMPI)
  set(ENV{CC}  cc) # C compiler for parallel build
  set(ENV{CXX} CC) # C++ compiler for parallel build
elseif()
  set(ENV{CC}  icc ) # C compiler for parallel build
  set(ENV{CXX} icpc) # C++ compiler for serial build
endif()

set(USER_CXX_FLAGS "-restrict -DMPICH_IGNORE_CXX_SEEK -std=c++14 -fopenmp")
set(USER_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -xHOST -fno-alias -restrict")
set(USER_CXX_FLAGS_DEBUG "-traceback -check=conversions,stack,uninit -check-pointers=rw -check-pointers-dangling=all-check-pointers-undimensioned -fp-stack-check -fp-trap=common -fp-trap-all=common")

set(FFTW_INCLUDE_DIR   $ENV{FFTW_INC})
set(FFTW_LIB           $ENV{FFTW_DIR}/libfftw3.a)
set(NETCDF_INCLUDE_DIR $ENV{NETCDF4_INCLUDE})
set(NETCDF_LIB_C       $ENV{NETCDF4_DIR}/lib/libnetcdf.a)
set(HDF5_LIB_1         $ENV{HDF5_DIR}/lib/libhdf5.a)
set(HDF5_LIB_2         $ENV{HDF5_DIR}/lib/libhdf5_hl.a)
set(SZIP_LIB           $ENV{SZIP_LIB}/libsz.a)
set(BOOST_INCLUDE_DIR  $ENV{BOOST_INCLUDE_DIR})

set(LIBS ${FFTW_LIB} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=restrict)
