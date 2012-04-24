FILE(GLOB sourcefiles "../src/*.cxx")
include_directories(${INCLUDE_DIRS} "../include")

add_executable(testmpi ${sourcefiles} testmpi.cxx mpicheck.cxx)
target_link_libraries(testmpi ${LIBS} m)
