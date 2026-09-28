cmake_minimum_required(VERSION 3.10)
project(PMTop)

set(CMAKE_CXX_COMPILER icpx)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost ")
set(CMAKE_CXX_FLAGS_DEBUG  "-O0 -DNDEBUG -g")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(MKL_LINK dynamic)
set(MKL_THREADING intel_thread)
set(MKL_INTERFACE lp64)
set(MKL_MPI intelmpi)

find_package(MKL REQUIRED)

if (MKL_FOUND)
    add_definitions(-D__MKL_MODULE)
endif()


find_package(arpackng REQUIRED)
if (arpackng_FOUND)
    message(STATUS "ARPACKNG found")
    add_definitions(-D__ARPACK_MODULE)
else()
    message(STATUS "ARPACKNG not found")
endif()

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found. Version: ${OpenMP_CXX_VERSION}")
    # 这一步通常在 target_link_libraries 中做，
    # 但如果你想全局开启，也可以在这里加 flags (不推荐，推荐在 src/CMakeLists.txt 里链接)
else()
    message(FATAL_ERROR "OpenMP not found! Parallelization will not work.")
endif()

find_package(nlohmann_json 3.11.2 REQUIRED)
if(nlohmann_json_FOUND)
    message(STATUS "nlohmann found")
    add_definitions(-D__NLOHMANN_MODULE)
endif()

find_package(spdlog REQUIRED)
if(spdlog_FOUND)
    message(STATUS "spdlog found")
    add_definitions(-D__SPDLOG_MODULE)
endif ()
add_subdirectory(src)

