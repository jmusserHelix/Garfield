cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(DelayedSignals)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# --- Build executables
add_executable(analytic analytic.C)
target_link_libraries(analytic Garfield::Garfield)

