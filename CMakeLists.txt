cmake_minimum_required(VERSION 3.18)

project(DFSweep) 

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 设置默认构建类型
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    add_definitions(-DENABLE_OMP)
    message(STATUS "OpenMP found, enabling OpenMP acceleration")
else()
    message(STATUS "OpenMP not found")
endif()


# 设置包含目录
file(GLOB INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd/MeshLib
    ${CMAKE_CURRENT_SOURCE_DIR}/cuh/
)

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

file(GLOB SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd/MeshLib/*.cpp
)


if(CUDA_AVAILABLE)
    file(GLOB CUDA_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/cu/*.cu
    )
    if(CUDA_SOURCE_FILES)
        list(APPEND SOURCE_FILES ${CUDA_SOURCE_FILES})
        set_source_files_properties(${CUDA_SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
        message(STATUS "CUDA source files: ${CUDA_SOURCE_FILES}")
    else()
        message(STATUS "No CUDA source files found")
    endif()
endif()


find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
    message("Eigen found")
    INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS
    "Eigen not found in system. Downloading Eigen using FetchContent...")
    file(DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/latest/download/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
    include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

    CPMAddPackage(
  NAME          Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG       3.4.0
)
    INCLUDE_DIRECTORIES(${Eigen_SOURCE_DIR})
endif()



message(STATUS "Downloading and building dependencies...")

include(FetchContent)

FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        v2.3.0
)
FetchContent_MakeAvailable(polyscope)

FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(CLI11)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} PRIVATE
	Eigen3::Eigen
    polyscope
    CLI11::CLI11)


if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()


install(TARGETS ${PROJECT_NAME}
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin
)
