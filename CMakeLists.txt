cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW) # Enforce the modern MSVC runtime policy
cmake_policy(SET CMP0174 NEW) # Suppress empty-keyword warning in vcpkg FindFFMPEG

project(N-Body C CXX)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE
  Main.cpp
  Simulation.cpp
  BHTree.cpp
  VideoRecorder.cpp
  glad.c
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

find_package(SDL3 CONFIG REQUIRED)
find_package(OpenGL REQUIRED)
find_package(FFMPEG REQUIRED)
target_include_directories(${PROJECT_NAME} PRIVATE ${FFMPEG_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE
  SDL3::SDL3
  OpenGL::GL
  ${FFMPEG_LIBRARIES}
)