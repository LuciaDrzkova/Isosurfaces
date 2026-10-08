cmake_minimum_required(VERSION 3.18)
project(IsosurfacesResearch LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)
include(CTest)

enable_testing()
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_VIS ON CACHE BOOL "" FORCE)
set(PMP_INSTALL ON CACHE BOOL "" FORCE)
set(PMP_STRICT_COMPILATION OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  pmp-library
  GIT_REPOSITORY https://github.com/pmp-library/pmp-library.git
  GIT_TAG 3.0.0
)
FetchContent_MakeAvailable(pmp-library)
if(TARGET imgui)
  file(GLOB _glfw_dirs "${pmp-library_SOURCE_DIR}/external/glfw-*")
  list(GET _glfw_dirs 0 _glfw_dir)
  target_include_directories(imgui PRIVATE "${_glfw_dir}/include")
endif()

add_library(iso_research INTERFACE)
target_include_directories(iso_research INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(iso_research INTERFACE pmp)

find_package(Threads REQUIRED)

add_executable(Isosurfaces src/main.cpp)
target_include_directories(Isosurfaces PRIVATE src include)
target_link_libraries(Isosurfaces PRIVATE iso_research pmp_vis Threads::Threads)

# Research-data executable:
# - compares the conforming LocalUnfolder baseline with
#   independently extracted local voxel MC patches,
# - stitches the voxel patch back into the global MC mesh,
# - sweeps c/resolution/radius/levels and records topology + scalar statistics,
# - supports the original cone, a non-quadric quartic saddle, and a sphere.
add_executable(IsosurfacesSweep src/research_sweep.cpp)
target_include_directories(IsosurfacesSweep PRIVATE include src)
target_link_libraries(IsosurfacesSweep PRIVATE iso_research pmp Threads::Threads)

# Single-configuration research viewer for the new voxel-patch -> global
# stitching path. The existing Isosurfaces executable remains the baseline GUI.
add_executable(IsosurfacesStitched src/stitched_viewer.cpp)
target_include_directories(IsosurfacesStitched PRIVATE include src)
target_link_libraries(IsosurfacesStitched PRIVATE iso_research pmp_vis Threads::Threads)

add_executable(IsosurfacesTests tests/tests.cpp)
target_include_directories(IsosurfacesTests PRIVATE include src)
target_link_libraries(IsosurfacesTests PRIVATE iso_research pmp Threads::Threads)

add_test(NAME IsosurfacesTests COMMAND IsosurfacesTests)
