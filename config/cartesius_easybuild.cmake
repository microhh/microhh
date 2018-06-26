# Cartesius using EasyBuild.
# Tested with:
#####
# GCC compiler (if USECUDA is on, build on GPU):
# module purge
# module load eb #(Use the new software development and installation framework EasyBuild currently implemented by SURFsara)
# module load surfsara
# module load CMake/3.9.5-GCCcore-6.4.0 #(Loads GCCcore as well)
# module load cuDNN/7.0.5-CUDA-9.0.176 #(Loads CUDA as well,cuDNN needed for Tensorflow-gpu)
# module load netCDF/4.5.0-foss-2017b #(Loads as well HDF5,cURL,sZIP,openMPI,FFTW3,GCC)
# module load netCDF-C++4/4.3.0-foss-2017b
# module load Doxygen/1.8.13-GCCcore-6.4.0
# module unload ScaLAPACK/2.0.2-gompi-2017b-OpenBLAS-0.2.20 #(Prevent crash during compiling: "(..)/microhh/src/../src/tools.cu([number]): (..) identifier [name] is undefined")
#####
# Intel Compiler (in all other cases):
# module purge
# module load eb (Use the new software development and installation framework EasyBuild, currently implemented by SURFsara)
# module load surfsara
# module load intel/2018a
# module load CMake/3.7.2-intel-2016b (Loads as well the MPI and Intel Compiler)
# module load cuDNN/7.0.5-CUDA-9.0.176 (Loads CUDA as well,cuDNN needed for Tensorflow-gpu)
# module load netCDF/4.4.1.1-intel-2016b (Loads as well HDF5,cURL,sZIP)
# module load netCDF-C++4/4.3.0-intel-2016b
# module load FFTW/3.3.5-intel-2016b
# module load Doxygen/1.8.11-intel-2016b

if(USEMPI)
    set(ENV{CC}  mpiicc ) # C compiler for parallel build
    set(ENV{CXX} mpiicpc) # C++ compiler for parallel build
    set(USEICC TRUE)
else()
    set(ENV{CC}  icc ) # C compiler for parallel build
    set(ENV{CXX} icpc) # C++ compiler for serial build
    set(USEICC TRUE)
endif()

if(USECUDA)
    set(ENV{CC}  gcc) # C compiler for serial build
    set(ENV{CXX} g++) # C++ compiler for serial build
    set(USEICC FALSE)
endif()

if(USEICC)
    set(USER_CXX_FLAGS_RELEASE "-Ofast -xAVX -axCORE-AVX-I,CORE-AVX2,CORE-AVX512")
else()
    set(USER_CXX_FLAGS_RELEASE "-Ofast -march=ivybridge") # -march optimized for the CPU present in Cartesius GPU nodes
endif()

set(USER_CXX_FLAGS "-std=c++14")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

#set(FFTW_INCLUDE_DIR "/home/bstratum/tools/fftw3_linked/include")
#set(FFTW_LIB         "/home/bstratum/tools/fftw3_linked/lib/libfftw3.a")
#set(HDF5_LIB_1 "/home/bstratum/tools/hdf5-1.8.17-gcc480/lib/libhdf5.a")
#set(HDF5_LIB_2 "/home/bstratum/tools/hdf5-1.8.17-gcc480/lib/libhdf5_hl.a")

set(FFTW_LIB       "fftw3")
set(FFTWF_LIB      "fftw3f")
set(NETCDF_LIB_C   "netcdf")
set(NETCDF_LIB_CPP "netcdf_c++4")
set(IRC_LIB        "irc")
set(HDF5_LIB       "hdf5")
set(SZIP_LIB       "sz")

set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB} ${SZIP_LIB} ${IRC_LIB} m z curl)
#set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR} ${NETCDF_INCLUDE_CXX_DIR})

if(USECUDA)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    #set(CUDA_LIB "cuda")
    #set(LIBS ${LIBS} ${CUDA_LIB})
    set(LIBS ${LIBS} -rdynamic $ENV{EBROOTCUDA}/lib64/libcufft.so)
    set(USER_CUDA_NVCC_FLAGS "-arch=sm_35")
    list(APPEND CUDA_NVCC_FLAGS "-std=c++14")
    list(APPEND CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
endif()

add_definitions(-DRESTRICTKEYWORD=__restrict__)
