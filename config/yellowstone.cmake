# Yellowstone
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

set(FFTW_INCLUDE_DIR   "/glade/apps/opt/fftw/3.3.4/intel/12.1.5/include")
set(FFTW_LIB           "/glade/apps/opt/fftw/3.3.4/intel/12.1.5/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/glade/apps/opt/netcdf/4.3.3.1/intel/12.1.5/include")
set(NETCDF_LIB_C       "/glade/apps/opt/netcdf/4.3.3.1/intel/12.1.5/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/glade/apps/opt/netcdf/4.3.3.1/intel/12.1.5/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/glade/apps/opt/hdf5/1.8.12/intel/12.1.5/lib/libhdf5.a")
set(HDF5_LIB_2         "/glade/apps/opt/hdf5/1.8.12/intel/12.1.5/lib/libhdf5_hl.a")
set(SZIP_LIB           "/glade/apps/opt/szip/2.1/intel/12.1.5/lib/libsz.a")
set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)

add_definitions(-DRESTRICTKEYWORD=restrict)
