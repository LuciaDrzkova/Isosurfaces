cmake_minimum_required(VERSION 3.16.3)
project(Isosurfaces LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # so VS Code's IntelliSense knows the flags/includes

include(FetchContent)

# ---------------------------------------------------------------------------
# pmp-library
# Brings its own vendored Eigen, GLFW, GLAD and ImGui, so nothing needs to be
# installed system-wide on either Windows or macOS besides a compiler and
# CMake itself. PMP_BUILD_VIEWERS builds the pmp::MeshViewer class we need.
# ---------------------------------------------------------------------------
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_VIS ON CACHE BOOL "" FORCE) # builds pmp_vis (MeshViewer, GLFW/ImGui window)
set(PMP_INSTALL ON CACHE BOOL "" FORCE) # NOTE: must stay ON — it's what makes pmp's own
                                         # headers visible to us at all (see CMakeLists.txt
                                         # comment below); it does not itself install anything
set(PMP_STRICT_COMPILATION OFF CACHE BOOL "" FORCE) # don't fail the build on pmp's own warnings

FetchContent_Declare(
  pmp-library
  GIT_REPOSITORY https://github.com/pmp-library/pmp-library.git
  GIT_TAG 3.0.0 # pinned: main has since renamed pmp/visualization/ to pmp/viewers/,
                # which would break MyViewer.h's #include paths
)
FetchContent_MakeAvailable(pmp-library)

# --- Workaround for a pmp-library 3.0.0 CMake bug -----------------------
# pmp-library's bundled imgui/CMakeLists.txt looks for GLFW's headers via
# ${CMAKE_SOURCE_DIR}/<path>, which is meant to be pmp-library's own source
# dir when it's built standalone. But CMAKE_SOURCE_DIR always resolves to
# the top-level project's dir (ours), so once pulled in via FetchContent
# it looks in the wrong place and fails with "GLFW/glfw3.h: No such file
# or directory". We fix it by adding the correct include path ourselves.
if(TARGET imgui)
  file(GLOB _glfw_dirs "${pmp-library_SOURCE_DIR}/external/glfw-*")
  list(GET _glfw_dirs 0 _glfw_dir)
  target_include_directories(imgui PRIVATE "${_glfw_dir}/include")
endif()
# --------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# autodiff
# Header-only (we only use the forward-mode "dual" API), so we fetch the
# source and just add it to the include path rather than building it.
# ---------------------------------------------------------------------------
FetchContent_Declare(
  autodiff
  GIT_REPOSITORY https://github.com/autodiff/autodiff.git
  GIT_TAG main
)
FetchContent_GetProperties(autodiff)
if(NOT autodiff_POPULATED)
  FetchContent_Populate(autodiff)
endif()

# ---------------------------------------------------------------------------
# Isosurfaces executable
# ---------------------------------------------------------------------------
add_executable(Isosurfaces src/main.cpp src/MyViewer.cpp)

target_include_directories(Isosurfaces PRIVATE
  ${CMAKE_SOURCE_DIR}/include
  ${autodiff_SOURCE_DIR}
)

# pmp     = core mesh data structure/algorithms
# pmp_vis = MeshViewer/TrackballViewer and the GLFW/ImGui window code
target_link_libraries(Isosurfaces PRIVATE pmp pmp_vis)

if(MSVC)
  target_compile_definitions(Isosurfaces PRIVATE _USE_MATH_DEFINES)
endif()

# Copy everything in assets/ (your .obj input meshes) next to the built
# executable after every build, so relative paths like "wavy_network.obj"
# in main.cpp resolve correctly no matter where CMake puts the binary.
add_custom_command(TARGET Isosurfaces POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_SOURCE_DIR}/assets
          $<TARGET_FILE_DIR:Isosurfaces>
)
