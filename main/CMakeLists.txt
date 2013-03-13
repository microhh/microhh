include_directories(${INCLUDE_DIRS} "../include")

add_executable(microhh microhh.cxx)
target_link_libraries(microhh microhhc ${LIBS} m)
