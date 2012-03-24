FILE(GLOB sourcefiles "*.cxx")
include_directories(${INCLUDE_DIRS})
add_executable(${PROJECT_NAME} ${sourcefiles})
target_link_libraries(${PROJECT_NAME} ${LIBS} m)
