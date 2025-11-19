cmake_minimum_required(VERSION 3.20)

project(DFSweep LANGUAGES CXX)

# 基本设置
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(MSVC)
    foreach(flag_var
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE)
        string(REPLACE "/MT" "/MD" ${flag_var} "${${flag_var}}")
    endforeach()
endif()

option(ENABLE_CUDA "Enable CUDA support" ON)
set(USE_CUDA FALSE)

if(WIN32 AND MSVC AND ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit QUIET)

    if(CUDAToolkit_FOUND)
        message(STATUS "CUDA toolkit found: ${CUDAToolkit_ROOT}")

        file(WRITE ${CMAKE_BINARY_DIR}/check_cuda_device.cu
        "#include <cuda_runtime.h>
        int main() {
            int count = 0;
            cudaError_t err = cudaGetDeviceCount(&count);
            return (err != cudaSuccess || count == 0) ? 1 : 0;
        }")

        try_compile(HAS_CUDA_DEVICE
            ${CMAKE_BINARY_DIR}
            ${CMAKE_BINARY_DIR}/check_cuda_device.cu
            CMAKE_FLAGS -DCMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}
            OUTPUT_VARIABLE CUDA_CHECK_OUTPUT
        )

        if(HAS_CUDA_DEVICE)
            message(STATUS "CUDA device detected, enabling CUDA support")
            add_definitions(-DENABLE_CUDA)
            set(USE_CUDA TRUE)

            project(DFSweep LANGUAGES CXX CUDA)
        else()
            message(WARNING "No CUDA device found. CUDA disabled.")
        endif()
    else()
        message(WARNING "CUDA toolkit not found. CUDA disabled.")
    endif()
else()
    message(STATUS "CUDA not checked (only active on Windows + MSVC).")
endif()

set(USE_OMP FALSE)

if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC|GNU|Clang")
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found")
        add_definitions(-DENABLE_OMP)
        set(USE_OMP TRUE)
    else()
        message(STATUS "OpenMP not found, continue without it")
    endif()
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd/MeshLib
)

if(USE_CUDA)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cuh)
endif()

file(GLOB SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd/MeshLib/*.cpp
)

if(USE_CUDA)
    file(GLOB CUDA_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/cu/*.cu
    )
    if(CUDA_SOURCE_FILES)
        list(APPEND SOURCE_FILES ${CUDA_SOURCE_FILES})
        set_source_files_properties(${CUDA_SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
        message(STATUS "CUDA sources added: ${CUDA_SOURCE_FILES}")
    endif()
endif()

find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS "Eigen not found, using FetchContent")
    include(FetchContent)

    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(Eigen)
    include_directories(${Eigen_SOURCE_DIR})
endif()

include(FetchContent)

FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG v2.3.0
)
FetchContent_MakeAvailable(polyscope)

FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(CLI11)

find_package(Doxygen QUIET)
find_package(Graphviz QUIET)
if(Doxygen_FOUND AND Graphviz_FOUND)
    message(STATUS "Doxygen found: ${DOXYGEN_EXECUTABLE}")
    message(STATUS "Graphviz found: ${Graphviz_DIR}")
    message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
    include(FetchContent)
    FetchContent_Declare(
    doxygen-awesome-css
    URL https://github.com/jothepro/doxygen-awesome-css/archive/refs/heads/main.zip
)
    FetchContent_MakeAvailable(doxygen-awesome-css)

    # Save the location the files were cloned into
    # This allows us to get the path to doxygen-awesome.css
    FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

    set(DOXYGEN_INPUT "${CMAKE_SOURCE_DIR}")
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/docs")
    set(DOXYGEN_EXTRACT_ALL YES)
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_GENERATE_LATEX YES)
    set(DOXYGEN_HAVE_DOT YES)
    set(DOXYGEN_CALL_GRAPH YES)
    set(DOXYGEN_REFERANCE_RELATION YES)
    set(DOXYFILE_CONTENT "
    PROJECT_NAME           = ${PROJECT_NAME}
    INPUT                  = ${DOXYGEN_INPUT} \
        ${CMAKE_SOURCE_DIR}/3rd/MeshLib \
        ${CMAKE_SOURCE_DIR}/include \
        ${CMAKE_SOURCE_DIR}/src \
        ${CMAKE_SOURCE_DIR}/cuh \
        ${CMAKE_SOURCE_DIR}/cu
    OUTPUT_DIRECTORY       = ${DOXYGEN_OUTPUT_DIRECTORY}
    EXTRACT_ALL            = ${DOXYGEN_EXTRACT_ALL}
    HAVE_DOT               = ${DOXYGEN_HAVE_DOT}
    GENERATE_TREEVIEW      = YES
    CALL_GRAPH             = ${DOXYGEN_CALL_GRAPH}
    REFERANCE_RELATION     = ${DOXYGEN_REFERANCE_RELATION}
    GENERATE_HTML          = ${DOXYGEN_GENERATE_HTML}
    GENERATE_LATEX         = ${DOXYGEN_GENERATE_LATEX}
    UML_LOOK               = ${Graphviz_FOUND}
    HTML_EXTRA_STYLESHEET  = ${AWESOME_CSS_DIR}/doxygen-awesome.css
    ")
    file(WRITE ${CMAKE_BINARY_DIR}/Doxyfile "${DOXYFILE_CONTENT}")
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating documentation with Doxygen"
    )
    install(DIRECTORY ${CMAKE_BINARY_DIR}/docs/html
            DESTINATION share/doc/${PROJECT_NAME}
    )
else()
    message(STATUS "Doxygen not found. Skipping documentation generation.")
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} PRIVATE
    Eigen3::Eigen
    polyscope
    CLI11::CLI11
)

if(USE_OMP)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS ${PROJECT_NAME}
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/docs/html
    DESTINATION share/doc/${PROJECT_NAME}
)


message(STATUS "CMake configuration finished.")
