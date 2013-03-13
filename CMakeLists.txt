INCLUDE(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
cmake_minimum_required (VERSION 2.8.1)
project (microhhc)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


#  make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE )
endif ()

#  make sure that the default system setting is macbook
if (NOT SYST)
  set (SYST macbook CACHE STRING
      "Choose the location: mpipc thunder fedora archlinux supermuc jugene" FORCE )
endif ()


# Set MPI on or off. By default it is on.
set (MPI TRUE CACHE BOOL "Choose to enable MPI or not.")
if(MPI)
 MESSAGE(STATUS "MPI is Enabled.")
 ADD_DEFINITIONS("-DPARALLEL")  # Set precompiler flag
 set (CXX_COMPILER_WRAPPER mpicxx)
 set (C_COMPILER_WRAPPER mpicc)
else()
  MESSAGE(STATUS "MPI is Disabled.")
endif()


get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER_ID} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER_ID} NAME)
if (C_COMPILER_NAME STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
elseif (C_COMPILER_NAME STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "-restrict")
  set(CMAKE_CXX_FLAGS_RELEASE "-xAVX -O3")
  set(CMAKE_CXX_FLAGS_DEBUG "-traceback -check=conversions,stack,uninit -check-pointers=rw -check-pointers-dangling=all-check-pointers-undimensioned -fp-stack-check -fp-trap=common -fp-trap-all=common ")
elseif (C_COMPILER_NAME STREQUAL "XL")
  if(MPI)
    set (C_COMPILER_WRAPPER mpcc)
    set (CXX_COMPILER_WRAPPER mpCC)
  endif()
  set(CMAKE_CXX_FLAGS_RELEASE "-qarch=pwr6 -qtune=pwr6 -O3 -qhot=simd -qenablevmx")
endif ()
# Set the compiler name to the wrapper, if necessary
MESSAGE(STATUS "Compiler " ${CMAKE_CXX_COMPILER})
MESSAGE(STATUS "Compiler " ${C_COMPILER_WRAPPER})


if (NOT ${C_COMPILER_WRAPPER}  STREQUAL "")
  CMAKE_FORCE_C_COMPILER(${C_COMPILER_WRAPPER} ${C_COMPILER_NAME})
  CMAKE_FORCE_CXX_COMPILER(${CXX_COMPILER_WRAPPER} ${CXX_COMPILER_NAME})
endif ()
MESSAGE(STATUS "Compiler " ${CMAKE_CXX_COMPILER})

if (SYST STREQUAL "mpipc")
  # MPI PC
  set(FFTW_INCLUDE_DIR   "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/include")
  set(FFTW_LIB           "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/lib/libfftw3.a")
  set(NETCDF_INCLUDE_DIR "/sw/squeeze-x64/netcdf-latest-static-gcc46/include")
  set(NETCDF_LIB_C       "/sw/squeeze-x64/netcdf-latest-static-gcc46/lib/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/sw/squeeze-x64/netcdf-latest-static-gcc46/lib/libnetcdf_c++.a")
  set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5.a")
  set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5_hl.a")
  set(SZIP_LIB           "/sw/squeeze-x64/szip-latest-static/lib/libsz.a")
  set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
elseif (SYST STREQUAL "blizzard")
  # Thunder
  set(FFTW_INCLUDE_DIR   "/pf/zmaw/m300041/local/include")
  set(FFTW_LIB           "/pf/zmaw/m300041/local/lib/libfftw3.a")
  set(NETCDF_INCLUDE_DIR "/sw/aix61/netcdf-4.1.3/include")
  set(NETCDF_LIB_C       "/sw/aix61/netcdf-4.1.3/lib/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/sw/aix61/netcdf-4.1.3/lib/libnetcdf_c++.a")
  set(HDF5_LIB_1         "/sw/aix61/hdf5-1.8.8/lib/libhdf5.a")
  set(HDF5_LIB_2         "/sw/aix61/hdf5-1.8.8/lib/libhdf5_hl.a")
  set(SZIP_LIB           "/sw/aix53/szip-2.1/lib/libsz.a")
  set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z)
elseif (SYST STREQUAL "thunder")
  # Thunder
  set(FFTW_INCLUDE_DIR   "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/include")
  set(FFTW_LIB           "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/lib/libfftw3.a")
  set(NETCDF_INCLUDE_DIR "/sw/squeeze-x64/netcdf-latest-static-gcc46/include")
  set(NETCDF_LIB_C       "/sw/squeeze-x64/netcdf-latest-static-gcc46/lib/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/sw/squeeze-x64/netcdf-latest-static-gcc46/lib/libnetcdf_c++.a")
  set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5.a")
  set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-latest-static/lib/libhdf5_hl.a")
  set(SZIP_LIB           "/sw/squeeze-x64/szip-latest-static/lib/libsz.a")
  set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
elseif (SYST STREQUAL "macbook")
  # Macbook
  set(FFTW_INCLUDE_DIR   "/opt/local/include")
  set(FFTW_LIB           "/opt/local/lib/libfftw3.a")
  set(NETCDF_INCLUDE_DIR "/opt/local/include")
  set(NETCDF_LIB_C       "/opt/local/lib/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/opt/local/lib/libnetcdf_c++.a")
  set(HDF5_LIB_1         "/opt/local/lib/libhdf5.a")
  set(HDF5_LIB_2         "/opt/local/lib/libhdf5_hl.a")
  set(SZIP_LIB           "")
  set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
elseif (SYST STREQUAL "fedora")
  # Fedora
  set(FFTW_INCLUDE_DIR   "/usr/include")
  set(FFTW_LIB           "/usr/lib64/libfftw3.a")
  set(NETCDF_INCLUDE_DIR "/usr/include")
  set(NETCDF_LIB_C       "/usr/lib64/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/usr/lib64/libnetcdf_c++.a")
  set(HDF5_LIB_1         "/usr/lib64/libhdf5.a")
  set(HDF5_LIB_2         "/usr/lib64/libhdf5_hl.a")
  set(SZIP_LIB           "")
  set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
elseif (SYST STREQUAL "supermuc")
  # SuperMUC
  set(FFTW_INCLUDE_DIR   "/lrz/sys/libraries/fftw/3.3.2/sse/include")
  set(FFTW_LIB           "/lrz/sys/libraries/fftw/3.3.2/sse/lib/libfftw3.a")
  set(NETCDF_INCLUDE_DIR "/lrz/sys/libraries/netcdf/4.2.1.1/include")
  set(NETCDF_LIB_C       "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf.a")
  set(NETCDF_LIB_CPP     "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf_c++.a")
  set(HDF5_LIB_1         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5.a")
  set(HDF5_LIB_2         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5_hl.a")
  set(SZIP_LIB           "/lrz/sys/libraries/hdf5/szip_2.1_u1/lib/libsz.a")
  set(LIBS ${FFTW_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
endif()

set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(main)
add_subdirectory(test)


###documentation
# add_custom_target(todo ALL)
# ADD_CUSTOM_COMMAND(TARGET todo POST_BUILD
#                   COMMAND echo "MICROHH TODO LIST" > TODO
#                   COMMAND date  >> TODO
#                   COMMAND grep -Rin \\todo  src | sed 's/!.*TODO//I' >>  TODO
#                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                   )
# find_package(Doxygen)
# IF (DOXYGEN_FOUND)
#   CONFIGURE_FILE(config/doxygen.conf.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
#   ADD_CUSTOM_TARGET(doc
#           ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
#           COMMENT "Generating API documentation with Doxygen" VERBATIM
#           DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
#   add_custom_target(install-doc
#           cp -r ${CMAKE_CURRENT_BINARY_DIR}/doc ${CMAKE_SOURCE_DIR}
#           DEPENDS doc)
# #   add_dependencies(doc-install doc)
# 
#   install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc DESTINATION ${CMAKE_SOURCE_DIR})
# ELSE (DOXYGEN)
#     MESSAGE(STATUS "WARNING: Doxygen not found - Reference manual will
# not be created")
# ENDIF ()
