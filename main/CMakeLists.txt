FILE(GLOB sourcefiles "../src/*.cxx")
include_directories(${INCLUDE_DIRS} "../include")

add_executable(init ${sourcefiles} init.cxx)
target_link_libraries(init ${LIBS} m)

add_executable(dns ${sourcefiles} dns.cxx)
target_link_libraries(dns ${LIBS} m)
