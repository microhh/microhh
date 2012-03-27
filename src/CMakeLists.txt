FILE(GLOB sourcefiles "*.cxx")
include_directories(${INCLUDE_DIRS} "../include")
add_executable(${PROJECT_NAME} ${sourcefiles})
target_link_libraries(${PROJECT_NAME} ${LIBS} m)
