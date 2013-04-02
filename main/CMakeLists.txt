include_directories(${INCLUDE_DIRS} "../include")

# retrieve the git hash from the current commit
exec_program(git ${CMAKE_CURRENT_SOURCE_DIR} ARGS "rev-parse --verify HEAD" OUTPUT_VARIABLE GITHASH ERROR_VARIABLE GITRETURN)
message(STATUS "Git hash " ${GITHASH})
message(STATUS "Git error" ${GITRETURN})

# send a precompiler statement replacing the git hash
add_definitions(-DGITHASH="${GITHASH}")
add_executable(microhh microhh.cxx)
target_link_libraries(microhh microhhc ${LIBS} m)
