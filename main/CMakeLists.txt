# 
#  MicroHH
#  Copyright (c) 2011-2023 Chiel van Heerwaarden
#  Copyright (c) 2011-2023 Thijs Heus
#  Copyright (c) 2014-2023 Bart van Stratum
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
  execute_process(
    COMMAND git describe --dirty
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE GITHASH
    RESULT_VARIABLE GITRETURN)
  if(NOT GITRETURN MATCHES 0)
    set(GITHASH "not available")
  endif()
else()
  set(GITHASH "not available")
endif()

message(STATUS "Git hash " ${GITHASH})

# send a precompiler statement replacing the git hash
add_definitions(-DGITHASH="${GITHASH}")

add_executable(microhh microhh.cxx)
target_include_directories(microhh PUBLIC ${INCLUDE_DIRS} "../include")

if(USECUDA)
    target_link_libraries(microhh microhhc curand rte_rrtmgp rte_rrtmgp_cuda rte_rrtmgp_cuda_rt rte_rrtmgp_kernels rte_rrtmgp_kernels_cuda rte_rrtmgp_kernels_cuda_rt ${LIBS} m)
else()
    target_link_libraries(microhh microhhc rte_rrtmgp_kernels rte_rrtmgp ${LIBS} m)
endif()
