cmake_minimum_required (VERSION 2.8.1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# load system specific settings if not set, force default.cmake
if(NOT SYST)
  set(SYST default)
endif()
include(${SYST} OPTIONAL RESULT_VARIABLE SYSTINC)

# trigger fatal error if illegal module is loaded
if(${SYSTINC} STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Config file config/" ${SYST} ".cmake does not exist.")
endif()

# set the default build type to RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
endif()

# start the project only after the system specific settings are loaded
project(microhhc CXX)

# enable or disable MPI, enabled by default
set(MPI TRUE CACHE BOOL "Choose to enable MPI or not.")
if(MPI)
  message(STATUS "MPI: enabled.")
  # add the necessary precompiler flag
  add_definitions("-DPARALLEL")
else()
  message(STATUS "MPI: disabled.")
endif()

# only set the compiler flags when the cache is created
# to enable editing of the flags in the CMakeCache.txt file
if(NOT HASCACHE)
  set(CMAKE_CXX_FLAGS ${USER_CXX_FLAGS} CACHE STRING
    "Flags used by the compiler during all build types." FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG ${USER_CXX_FLAGS_DEBUG} CACHE STRING
    "Flags used by the compiler during debug builds." FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE ${USER_CXX_FLAGS_RELEASE} CACHE STRING
    "Flags used by the compiler during release builds." FORCE)
  message(STATUS "Build Type: " ${CMAKE_BUILD_TYPE})
  set(HASCACHE TRUE CACHE BOOL "CMakeCache.txt created." FORCE)
  # make sure that ccmake only contains build type and MPI setting
  mark_as_advanced(HASCACHE)
  mark_as_advanced(CMAKE_INSTALL_PREFIX)
endif()

# print the compiler flags to the screen
if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  message(STATUS "Compiler flags: " ${CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS_RELEASE})
else()
  message(STATUS "Compiler flags: " ${CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS_DEBUG})
endif()

# set the include dirs
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(main)
add_subdirectory(test)

# documentation
add_custom_target(todo ALL)
# if a non-default sed is not set revert to standard sed
if(NOT DEFINED GNU_SED)
  set(GNU_SED "sed")
endif()
add_custom_command(TARGET todo POST_BUILD
  COMMAND echo "MICROHH TODO LIST" > TODO
  COMMAND date >> TODO
  COMMAND grep -Rin \\todo  src include main| ${GNU_SED} 's/ *\\/\\/.*TODO */ /I' >> TODO
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(config/doxygen.conf.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
    COMMENT "Generating API documentation with Doxygen." VERBATIM
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc DESTINATION ${CMAKE_SOURCE_DIR})
else(DOXYGEN)
  message(STATUS "WARNING: Doxygen not found - Reference manual will not be created.")
endif()
