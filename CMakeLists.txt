cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW) # Enforce the modern MSVC runtime policy

project(N-Body CXX)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE
  Main.cpp
  Simulation.cpp
  BHTree.cpp
  VectorMath.cpp
)

if(MSVC)
  #statically link MSVC runtime library
  set_target_properties(${PROJECT_NAME} PROPERTIES
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

  #set compiler optimization flags
  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Release>:/O2>
    $<$<CONFIG:Debug>:/Od>
  )
else() # For compilers like GCC or Clang
  #set compiler optimization flags
  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Release>:-O2>
    $<$<CONFIG:Debug>:-O0>
  )
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

find_library(OPENGLl32_LIB_FILE
  NAMES opengl32
)

if(OPENGLl32_LIB_FILE)
  target_link_libraries(${PROJECT_NAME} PRIVATE
    "${OPENGLl32_LIB_FILE}"
  )
else()
  message(FATAL_ERROR "Custom library not found. Please check the path.")
endif()