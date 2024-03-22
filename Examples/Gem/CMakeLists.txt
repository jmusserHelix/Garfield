# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(gem)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executable------------------------------------------------------------
add_executable(gem gem.C)
target_link_libraries(gem Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file ELIST.lis MPLIST.lis NLIST.lis PRNSOL.lis)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
