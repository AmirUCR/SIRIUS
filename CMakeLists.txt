cmake_minimum_required(VERSION 3.20)
project(SIRIUS)
set(CMAKE_CXX_STANDARD 17)

set(OR_TOOLS_DIR "${CMAKE_SOURCE_DIR}/ortools")
include_directories(${OR_TOOLS_DIR}/include)
link_directories(${OR_TOOLS_DIR}/lib)

add_executable(sirius
    sirius.cpp
    sirius/src/utils.cpp
    sirius/src/tables.cpp
    sirius/src/config.cpp
    sirius/src/solver.cpp
    sirius/src/logging.cpp
    sirius/src/instance.cpp
    sirius/src/interface.cpp
)

find_package(Threads REQUIRED)

target_include_directories(sirius PRIVATE
    "${CMAKE_SOURCE_DIR}/sirius/include"
)

target_link_libraries(sirius PRIVATE ortools Threads::Threads)

set_target_properties(sirius PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/sirius/build"
)
