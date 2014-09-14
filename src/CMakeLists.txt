# 
#  MicroHH
#  Copyright (c) 2011-2014 Chiel van Heerwaarden
#  Copyright (c) 2011-2014 Thijs Heus
#  Copyright (c)      2014 Bart van Stratum
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
FILE(GLOB sourcefiles "../src/*.cxx" "../src/*.cu")
include_directories(${INCLUDE_DIRS} "../include")

if(USECUDA)
  cuda_add_library(microhhc STATIC ${sourcefiles})
else()
  add_library(microhhc STATIC ${sourcefiles})
endif()
