cmake_minimum_required (VERSION 2.6)
project (microhhc)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(FFTW_INCLUDE_DIR "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/include")
set(FFTW_LIB         "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/lib/libfftw3.a")

set(FFTW_INCLUDE_DIR "/opt/local/include")
set(FFTW_LIB         "/opt/local/lib/libfftw3.dylib")

#set(FFTW_INCLUDE_DIR "/lrz/sys/libraries/fftw/3.3/include")
#set(FFTW_LIB         "/lrz/sys/libraries/fftw/3.3/lib/libfftw3.a")

#set(FFTW_INCLUDE_DIR "/bgsys/local/fftw3/3.3/fftw_450/include")
#set(FFTW_LIB         "/bgsys/local/fftw3/3.3/fftw_450/lib/libfftw3.a")

#set(FFTW_INCLUDE_DIR "/bgsys/local/fftw3/3.3.2/fftw_g/include")
#set(FFTW_LIB         "/bgsys/local/fftw3/3.3.2/fftw_g/lib/libfftw3.a")

set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

add_subdirectory(src)
add_subdirectory(main)
add_subdirectory(test)

