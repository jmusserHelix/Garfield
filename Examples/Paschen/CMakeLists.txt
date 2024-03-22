# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Paschen)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executables-----------------------------------------------------------
add_executable(gastable gastable.C)
target_link_libraries(gastable Garfield::Garfield)

add_executable(paschenariso paschenariso.C)
target_link_libraries(paschenariso Garfield::Garfield)
