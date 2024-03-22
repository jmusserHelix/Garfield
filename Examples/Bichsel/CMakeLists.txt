cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Bichsel)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executable------------------------------------------------------------
add_executable(test test.C)
target_link_libraries(test Garfield::Garfield)
