# Enables all the relevant warnings for the different compilers
function(garfield_enable_default_warnings target_name)
  target_compile_options(
    ${target_name}
    PRIVATE
      $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
      -Wall
      -Wextra
      -pedantic
      -Wfatal-errors
      -Wno-missing-braces
      -Wshadow>
      $<$<CXX_COMPILER_ID:MSVC>:/W3>)
endfunction()

# ########################################################################################
# Sets the typical compile options for a given target, besides the warnings which are set
# previously
function(garfield_enable_default_compiler_options target_name)
  target_compile_options(${target_name} PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/EHsc >)
endfunction()

# ########################################################################################
# There are two ways for a library to find its dependencies. They can be hardcoded in the
# library header or they can be relative paths the runtime linker should look for. More
# info here: https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling The
# following function allows to build relocatable libraries which do not encode the RPATH
# in the installed library so that it relies on the runtime linker for the linking

function(garfield_set_default_rpath_properties target_name)
  set_property(TARGET ${target_name} PROPERTY SKIP_BUILD_RPATH FALSE)
  set_property(TARGET ${target_name} PROPERTY BUILD_WITH_INSTALL_RPATH FALSE)
  set_property(TARGET ${target_name} PROPERTY INSTALL_RPATH "")
  set_property(TARGET ${target_name} PROPERTY INSTALL_RPATH_USE_LINK_PATH FALSE)
  set_property(TARGET ${target_name} PROPERTY MACOSX_RPATH TRUE)
endfunction()

# ########################################################################################
# The function sets all the default properties for the specified target
function(garfield_set_all_default_properties target_name)
  garfield_enable_default_warnings(${target_name})
  garfield_enable_default_compiler_options(${target_name})
  garfield_set_default_rpath_properties(${target_name})
endfunction()

# ########################################################################################
# Forces the color output for GNU and Clang in all cases Useful to enable that for Ninja
# build
function(force_color_output)
  option(force_colored_output "Always produce ANSI-colored output (GNU/Clang only)." TRUE)
  mark_as_advanced(force_colored_output)
  if(force_colored_output)
    add_compile_options(
      $<$<CXX_COMPILER_ID:GNU>:-fdiagnostics-color=always>
      $<$<AND:$<CXX_COMPILER_ID:Clang>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,4.0.0>>:-fcolor-diagnostics>
    )
  endif()
endfunction()
