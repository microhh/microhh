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

# Select correct compilers for GCC + parallel/serial:
if(USEMPI)
    set(ENV{CXX} mpicxx)
    set(ENV{FC} mpif90)
else()
    set(ENV{CXX} CC)
    set(ENV{FC} ftn)
endif()

# Set compiler flags / options.

#set(USER_CXX_FLAGS_RELEASE "-O3 -march=znver2 -mtune=znver2 -mfma -mavx2 -m3dnow -fomit-frame-pointer")
#set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
#
#set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
#set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=znver2 -mtune=znver2 -mfma -mavx2 -m3dnow -fomit-frame-pointer")

set(USER_CXX_FLAGS "")
set(USER_CXX_FLAGS_RELEASE "-O3")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(USER_FC_FLAGS "")
set(USER_FC_FLAGS_RELEASE "-O3")
set(USER_FC_FLAGS_DEBUG "-O0 -g -Wall")

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
