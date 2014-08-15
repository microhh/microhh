FILE(GLOB sourcefiles "../src/*.cxx" "../src/*.cu")
include_directories(${INCLUDE_DIRS} "../include")

if(USECUDA)
  cuda_add_library(microhhc STATIC ${sourcefiles})
else()
  add_library(microhhc STATIC ${sourcefiles})
endif()
