# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Voxel)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Define executables----------------------------------------------------------
add_executable(plotfield plotfield.C)
target_link_libraries(plotfield Garfield::Garfield)

add_executable(savefield savefield.C)
target_link_libraries(savefield Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file Efield.txt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
