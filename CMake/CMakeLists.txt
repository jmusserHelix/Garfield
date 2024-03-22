# Configuring the package for exporting
include(CMakePackageConfigHelpers)

# This is a valid destination for both Win and Unix
set(CMAKE_INSTALL_DIR share/${PROJECT_NAME}/cmake)

# Creating the target export file
install(
  EXPORT "${PROJECT_NAME}Targets"
  FILE "${PROJECT_NAME}Targets.cmake"
  NAMESPACE Garfield::
  DESTINATION "${CMAKE_INSTALL_DIR}")

# Creates the versioning file
set(VersionFilePath "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
write_basic_package_version_file(
  ${VersionFilePath}
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

# Prepares the configuration file
set(ConfigFileName "${PROJECT_NAME}Config.cmake")
set(ConfigFileTemplate ${CMAKE_CURRENT_LIST_DIR}/${ConfigFileName}.in)
configure_package_config_file(
  "${ConfigFileTemplate}" "${CMAKE_CURRENT_BINARY_DIR}/${ConfigFileName}"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DIR}")

# Installs the configuration file
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${ConfigFileName}" "${VersionFilePath}"
        DESTINATION "${CMAKE_INSTALL_DIR}")

find_package(Python)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in
               ${PROJECT_BINARY_DIR}/setupGarfield.sh @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.csh.in
               ${PROJECT_BINARY_DIR}/setupGarfield.csh @ONLY)

install(FILES ${PROJECT_BINARY_DIR}/setupGarfield.sh 
              ${PROJECT_BINARY_DIR}/setupGarfield.csh
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Garfield)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Garfield)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in
               ${PROJECT_BINARY_DIR}/Garfield/__init__.py @ONLY)
install(DIRECTORY ${PROJECT_BINARY_DIR}/Garfield/ 
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/Garfield)
