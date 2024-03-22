cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Electroluminescence)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Define executables---------------------------------------------------------
add_executable(example example.C)
target_link_libraries(example Garfield::Garfield)

add_executable(userHandle userHandle.C)
target_link_libraries(userHandle Garfield::Garfield)
