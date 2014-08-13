FILE(GLOB sourcefiles "../src/*.cxx" "../src/*.cu")
include_directories(${INCLUDE_DIRS} "../include")

cuda_add_library(microhhc STATIC ${sourcefiles})
