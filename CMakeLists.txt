set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)

INCLUDE(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
cmake_minimum_required (VERSION 2.8.1)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


#  make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE )
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
      "Choose the type of build, options are: None Debug Release."  FORCE)
   MESSAGE(STATUS${CMAKE_BUILD_TYPE} ${TEMP})

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

#Load system specific settings
if (NOT SYST)
  set (SYST default CACHE STRING
      "Choose the location: mpipc thunder fedora archlinux supermuc jugene" FORCE )
endif ()
include (${SYST} OPTIONAL)

#Start the project only after all the variables are set
project (microhhc)

get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER_ID} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER_ID} NAME)
MESSAGE(STATUS "Compiler " ${CMAKE_CXX_COMPILER})
MESSAGE(STATUS "Compiler " ${C_COMPILER_WRAPPER})


if (NOT ${C_COMPILER_WRAPPER}  STREQUAL "")
  CMAKE_FORCE_C_COMPILER(${C_COMPILER_WRAPPER} ${C_COMPILER_NAME})
  CMAKE_FORCE_CXX_COMPILER(${CXX_COMPILER_WRAPPER} ${CXX_COMPILER_NAME})
endif ()
MESSAGE(STATUS "Compiler " ${CMAKE_CXX_COMPILER})

set(CMAKE_CXX_FLAGS ${USER_CXX_FLAGS} )
set(CMAKE_CXX_FLAGS_RELEASE ${USER_CXX_FLAGS_RELEASE} )
set(CMAKE_CXX_FLAGS_DEBUG ${USER_CXX_FLAGS_DEBUG})
MESSAGE(STATUS "Build Type " ${CMAKE_BUILD_TYPE})
if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})
else()
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG})
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
