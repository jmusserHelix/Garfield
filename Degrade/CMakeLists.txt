set(libname "degrade")
set(lib_fullver 3.17)
set(lib_apiver 3)

# Creating the target as a shared library
message(STATUS "Creating a new library named ${libname} with API version ${lib_apiver}")
add_library(${libname} SHARED "")
add_library(Garfield::${libname} ALIAS ${libname})
target_sources(${libname} PRIVATE degrad.f)
target_compile_options(${libname} PRIVATE -std=legacy)
# Sets the API version of the library
set_property(TARGET ${libname} PROPERTY SOVERSION ${lib_apiver})

install(
  TARGETS ${libname}
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
