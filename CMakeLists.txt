cmake_minimum_required (VERSION 2.6)
project (microhhc)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_subdirectory(src)

set(FFTW_INCLUDE_DIR "/sw/lenny-x64/numerics/fftw-3.3-openmp-gccsys/include")
set(FFTW_LIB_DIR     "/sw/lenny-x64/numerics/fftw-3.3-openmp-gccsys/lib")

set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIB_DIRS     ${FFTW_LIB_DIR})

