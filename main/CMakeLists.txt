# 
#  MicroHH
#  Copyright (c) 2011-2017 Chiel van Heerwaarden
#  Copyright (c) 2011-2017 Thijs Heus
#  Copyright (c) 2014-2017 Bart van Stratum
# 
#  This file is part of MicroHH
# 
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#
include_directories(${INCLUDE_DIRS} "../include")

# retrieve the git hash from the current commit
find_package(Git)
if(GIT_EXECUTABLE)
  exec_program(git ${CMAKE_CURRENT_SOURCE_DIR} ARGS "describe --dirty" OUTPUT_VARIABLE GITHASH RETURN_VALUE GITRETURN)
  if(NOT GITRETURN MATCHES 0)
    set(GITHASH "not available")
  endif()
else()
  set(GITHASH "not available")
endif()

message(STATUS "Git hash " ${GITHASH})

# send a precompiler statement replacing the git hash
add_definitions(-DGITHASH="${GITHASH}")

# double precision version
if(USECUDA)
  cuda_add_executable(microhh microhh.cxx)
  target_link_libraries(microhh microhhc rrtmg ${LIBS} m)
else()
  add_executable(microhh microhh.cxx)
  target_link_libraries(microhh microhhc rrtmg ${LIBS} m)
endif()

# single precision version
if(USECUDA)
  cuda_add_executable(microhh_single microhh_single.cxx)
  target_link_libraries(microhh_single microhhc rrtmg ${LIBS} m)
else()
  add_executable(microhh_single microhh_single.cxx)
  target_link_libraries(microhh_single microhhc rrtmg ${LIBS} m)
endif()
