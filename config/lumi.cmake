
# *  ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒  *   *      *
#                                                       *      *  *
#   * ████       ████   ████   █████▄    ▄█████   ████     *     *
# *   ████       ████   ████   ████ █▄  ▄█ ████   ████         ,    *,
#     ████       ████   ████   ████  ████  ████   ████  *   *  |\_ _/|
#     ████       ████   ████   ████   ▀▀   ████   ████   *    .| ." ,|
#  *  ████       ████   ████   ████        ████   ████        /(  \_\)
#     ████       ████   ████   ████        ████   ████       /    ,-,|
# *   ████▄▄▄▄▄  ▀███   ███▀   ████        ████   ████ *    * /      \
#     █████████    ▀▀███▀▀     ████        ████   ████  * ,/  (      *
# *                                                     ,/       |  /
#  * ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒/    \  * || |
#                 *              *               ,_   (       )| || |
#*   *    *    The Supercomputer of the North  * | `\_|   __ /_| || |
#        **               *            * *       \_____\______)\__)__)
#   .********----------*******-------******----------****************.

# Programming environment: "cray" or "gnu"
set(PRGENV "cray")

# Select correct compilers based on PRGENV and parallel/serial:
if(USEMPI)
    set(ENV{CXX} mpicxx)
    set(ENV{FC} mpif90)
else()
    set(ENV{CXX} CC)
    set(ENV{FC} ftn)
endif()

# Set compiler flags / options based on PRGENV.
if(PRGENV STREQUAL "gnu")
    set(USER_CXX_FLAGS "-std=c++20")
    set(USER_CXX_FLAGS_RELEASE "-O3 -march=znver3 -mtune=znver3 -mfma -mavx2 -fomit-frame-pointer")
    set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

    set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
    set(USER_FC_FLAGS_RELEASE "-O3 -march=znver3 -mtune=znver3 -mfma -mavx2 -fomit-frame-pointer")
    set(USER_FC_FLAGS_DEBUG "-O0 -g")
else()
    # Cray compiler flags
    set(USER_CXX_FLAGS "-std=c++20")
    set(USER_CXX_FLAGS_RELEASE "-Ofast -ffp=3 -flto -funroll-loops")
    set(USER_CXX_FLAGS_DEBUG "-O0 -g")

    set(USER_FC_FLAGS "")
    set(USER_FC_FLAGS_RELEASE "-O3 -hfp3")
    set(USER_FC_FLAGS_DEBUG "-O0 -g")
endif()

# Set libraries.
link_directories(
    $ENV{FFTW_DIR}
    $ENV{NETCDF_DIR}/lib
    $ENV{HDF5_DIR}/lib
)

include_directories(
    $ENV{NETCDF_DIR}/include
    $ENV{FFTW_INC})

set(NETCDF_LIB_C "netcdf")
set(FFTW_LIB "fftw3")
set(FFTWF_LIB "fftw3f")
set(IRC_LIB "irc")
set(IRC_LIB "")
set(HDF5_LIB "hdf5")
set(SZIP_LIB "sz")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB} ${SZIP_LIB} ${IRC_LIB} m z curl)

add_definitions(-DRESTRICTKEYWORD=__restrict__)
add_definitions(-DDISABLE_2D_MPIIO=1)
add_definitions(-DRTE_USE_CBOOL)
