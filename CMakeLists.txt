cmake_minimum_required (VERSION 2.6)
project (microhhc)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(FFTW_INCLUDE_DIR "/sw/lenny-x64/numerics/fftw-3.3-openmp-gccsys/include")
set(FFTW_LIB         "/sw/lenny-x64/numerics/fftw-3.3-openmp-gccsys/lib/libfftw3.a")

#set(FFTW_INCLUDE_DIR "/opt/local/include")
#set(FFTW_LIB         "/opt/local/lib/libfftw3.dylib")

set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

add_subdirectory(src)
