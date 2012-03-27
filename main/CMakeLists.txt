FILE(GLOB sourcefiles "../src/*.cxx")
include_directories(${INCLUDE_DIRS} "../include")

add_executable(dns ${sourcefiles} dns.cxx)
target_link_libraries(dns ${LIBS} m)
