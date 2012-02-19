#include_directories(${MICROHH_INCLUDE_DIRS})
FILE(GLOB sourcefiles "*.cxx")
add_executable(${PROJECT_NAME} ${sourcefiles})
#target_link_libraries(${PROJECT_NAME} ${MICROHH_LIBRARIES})
