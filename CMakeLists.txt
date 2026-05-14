cmake_minimum_required(VERSION 3.24)

project(mold-check LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(VCLIB_BUILD_MODULE_EXTERNAL ON CACHE BOOL "Build VCLib external module")
set(VCLIB_BUILD_MODULE_RENDER ON CACHE BOOL "Build VCLib render module")

add_subdirectory(vclib)

add_executable(mold-check main.cpp)

target_link_libraries(mold-check PRIVATE
    vclib::external
    vclib::render)

target_compile_definitions(mold-check PRIVATE
    NOMINMAX
    _USE_MATH_DEFINES
    MESHES_PATH="${CMAKE_CURRENT_SOURCE_DIR}/meshes"
    RESULTS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/output")

add_custom_target(run
    COMMAND $<TARGET_FILE:mold-check>
    DEPENDS mold-check
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    USES_TERMINAL)
