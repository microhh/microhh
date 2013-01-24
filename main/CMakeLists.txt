include_directories(${INCLUDE_DIRS} "../include")

add_executable(init init.cxx)
target_link_libraries(init microhhc ${LIBS} m)

add_executable(dns dns.cxx)
target_link_libraries(dns microhhc ${LIBS} m)

add_executable(stats stats.cxx)
target_link_libraries(stats microhhc ${LIBS} m)
