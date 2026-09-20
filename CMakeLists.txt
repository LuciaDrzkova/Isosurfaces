cmake_minimum_required(VERSION 3.18)
project(Isosurfaces LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # so VS Code's IntelliSense knows the flags/includes

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)

# ---------------------------------------------------------------------------
# pmp-library
# Brings its own vendored Eigen, GLFW, GLAD and ImGui, so nothing needs to be
# installed system-wide besides a compiler and CMake. PMP_BUILD_VIS builds
# pmp::MeshViewer (the GLFW/ImGui window).
# ---------------------------------------------------------------------------
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_VIS ON CACHE BOOL "" FORCE)
set(PMP_INSTALL ON CACHE BOOL "" FORCE) # must stay ON: pmp only exposes its own
                                        # include directories when this is set
set(PMP_STRICT_COMPILATION OFF CACHE BOOL "" FORCE) # don't fail on pmp's own warnings

FetchContent_Declare(
  pmp-library
  GIT_REPOSITORY https://github.com/pmp-library/pmp-library.git
  GIT_TAG 3.0.0 # pinned: main renamed pmp/visualization/ to pmp/viewers/
)
FetchContent_MakeAvailable(pmp-library)

# Workaround for a pmp-library 3.0.0 bug: its bundled imgui/CMakeLists.txt
# looks for GLFW's headers under ${CMAKE_SOURCE_DIR}, which is *our* source
# dir when pmp is pulled in via FetchContent, so the build fails with
# "GLFW/glfw3.h: No such file or directory". Add the correct path ourselves.
if(TARGET imgui)
  file(GLOB _glfw_dirs "${pmp-library_SOURCE_DIR}/external/glfw-*")
  list(GET _glfw_dirs 0 _glfw_dir)
  target_include_directories(imgui PRIVATE "${_glfw_dir}/include")
endif()

# ---------------------------------------------------------------------------
# autodiff (header-only; only forward-mode dual numbers are used)
# Its own CMakeLists insists on a system-wide Eigen, which we don't need, so
# SOURCE_SUBDIR points at a directory that doesn't exist: MakeAvailable then
# only downloads the sources and we expose the headers ourselves.
# ---------------------------------------------------------------------------
FetchContent_Declare(
  autodiff
  GIT_REPOSITORY https://github.com/autodiff/autodiff.git
  GIT_TAG v1.1.2
  SOURCE_SUBDIR do-not-configure
)
FetchContent_MakeAvailable(autodiff)

add_library(autodiff INTERFACE)
add_library(autodiff::autodiff ALIAS autodiff)
target_include_directories(autodiff SYSTEM INTERFACE ${autodiff_SOURCE_DIR})

# ---------------------------------------------------------------------------
# The program. New .cpp files in src/method/ or src/app/ are picked up
# automatically (re-run "cmake -S . -B build", or just build again).
#
#   src/method/  the maths: surfaces, remeshing, measurements, research code
#   src/app/     the window, file browser, command line, saved settings
# ---------------------------------------------------------------------------
file(GLOB METHOD_SOURCES CONFIGURE_DEPENDS src/method/*.cpp)
file(GLOB APP_SOURCES CONFIGURE_DEPENDS src/app/*.cpp)

# "method" has no window code, so the tests can use it without OpenGL
add_library(isosurfaces_method STATIC ${METHOD_SOURCES})
target_include_directories(isosurfaces_method PUBLIC src)
target_link_libraries(isosurfaces_method PUBLIC pmp autodiff::autodiff)

if(MSVC)
  target_compile_options(isosurfaces_method PRIVATE /W4)
else()
  target_compile_options(isosurfaces_method PRIVATE -Wall -Wextra)
endif()

find_package(Threads REQUIRED) # the window runs jobs on a worker thread
add_executable(Isosurfaces ${APP_SOURCES})
target_link_libraries(Isosurfaces PRIVATE isosurfaces_method pmp_vis Threads::Threads)

# ---------------------------------------------------------------------------
# Tests (run with ctest)
# ---------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
  add_executable(isosurfaces_tests tests/tests.cpp src/app/Settings.cpp)
  target_link_libraries(isosurfaces_tests PRIVATE isosurfaces_method)
  add_test(NAME isosurfaces_tests COMMAND isosurfaces_tests)
endif()
