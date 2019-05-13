FILE(GLOB sourcefiles "../src_rrtmgp/*.cpp")
include_directories("../include_rrtmgp" SYSTEM ${INCLUDE_DIRS})

if(USECUDA)
  cuda_add_library(rrtmgp STATIC ${sourcefiles})
else()
  add_library(rrtmgp STATIC ${sourcefiles})
endif()
