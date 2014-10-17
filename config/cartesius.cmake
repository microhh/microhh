# Cartesius SurfSARA
if(USEMPI)
  set(ENV{CXX} mpiicpc) # compiler for parallel build
elseif()
  set(ENV{CXX} icpc) # compiler for serial build
endif()

set(USER_CXX_FLAGS "-restrict -DMPICH_IGNORE_CXX_SEEK")
set(USER_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -xHOST -fno-alias -restrict -vec-report1 -no-prec-div")
set(USER_CXX_FLAGS_DEBUG "-traceback -check=conversions,stack,uninit -check-pointers=rw -check-pointers-dangling=all-check-pointers-undimensioned -fp-stack-check -fp-trap=common -fp-trap-all=common")

set(FFTW_INCLUDE_DIR   "/hpc/sw/fftw3avx-3.3.3-intel-impi/include")
set(FFTW_LIB           "/hpc/sw/fftw3avx-3.3.3-intel-impi/lib/libfftw3.a")
set(NETCDF_INCLUDE_DIR "/hpc/sw/netcdf-4.1.3-intel-seq/include")
set(NETCDF_LIB_C       "/hpc/sw/netcdf-4.1.3-intel-seq/lib/libnetcdf.a")
set(NETCDF_LIB_CPP     "/hpc/sw/netcdf-4.1.3-intel-seq/lib/libnetcdf_c++.a")
set(HDF5_LIB_1         "/hpc/sw/hdf5-1.8.12-intel-seq/lib/libhdf5.a")
set(HDF5_LIB_2         "/hpc/sw/hdf5-1.8.12-intel-seq/lib/libhdf5_hl.a")
set(SZIP_LIB           "/hpc/sw/szip-2.1-intel//lib/libsz.a")

set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)

add_definitions(-DRESTRICTKEYWORD=restrict)
