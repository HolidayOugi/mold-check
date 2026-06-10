cmake_minimum_required(VERSION 3.24)

project(mold-check LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(VCLIB_BUILD_MODULE_EXTERNAL ON CACHE BOOL "Build VCLib external module")
set(VCLIB_BUILD_MODULE_RENDER ON CACHE BOOL "Build VCLib render module")

set(VCLIB_ALLOW_SYSTEM_CGAL OFF CACHE BOOL "Use system CGAL")
set(VCLIB_ALLOW_DOWNLOAD_CGAL ON CACHE BOOL "Download CGAL")

add_subdirectory(vclib)

add_executable(mold-check main.cpp)

target_link_libraries(mold-check PRIVATE
    vclib::external
    vclib::render
)

target_compile_definitions(mold-check PRIVATE
    NOMINMAX
    _USE_MATH_DEFINES
    MESHES_PATH="${CMAKE_CURRENT_SOURCE_DIR}/meshes"
    RESULTS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/output"
)

if(WIN32)
    set(VCLIB_CGAL_AUX_DLL_DIR
        "${CMAKE_BINARY_DIR}/_deps/cgal-src/auxiliary/gmp/bin"
    )

    set(EMBREE_DLL_DIR
        "C:/ProgramData/chocolatey/lib/embree.install/tools/bin"
    )

    if(EXISTS "${VCLIB_CGAL_AUX_DLL_DIR}")
        add_custom_command(TARGET mold-check POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E echo "Copying CGAL/GMP/MPFR DLLs..."
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${VCLIB_CGAL_AUX_DLL_DIR}"
                "$<TARGET_FILE_DIR:mold-check>"
        )
    else()
        message(WARNING "CGAL auxiliary DLL dir not found: ${VCLIB_CGAL_AUX_DLL_DIR}")
    endif()

    if(EXISTS "${EMBREE_DLL_DIR}")
        add_custom_command(TARGET mold-check POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E echo "Copying Embree DLLs..."
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${EMBREE_DLL_DIR}"
                "$<TARGET_FILE_DIR:mold-check>"
        )
    else()
        message(WARNING "Embree DLL dir not found: ${EMBREE_DLL_DIR}")
    endif()
endif()

add_custom_target(run
    COMMAND $<TARGET_FILE:mold-check>
    DEPENDS mold-check
    WORKING_DIRECTORY $<TARGET_FILE_DIR:mold-check>
    USES_TERMINAL
)