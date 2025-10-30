cmake_minimum_required(VERSION 3.18)  # keep compatible with manylinux images
project(SIRIUS VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Send all targets' outputs into the binary dir so setup.py can find them
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Os -ffunction-sections -fdata-sections)
  add_link_options(-Wl,--gc-sections -s)
endif()

include(FetchContent)
find_package(Threads REQUIRED)

# OR-Tools
FetchContent_Declare(
  or-tools
  GIT_REPOSITORY https://github.com/google/or-tools.git
  GIT_TAG        v9.11
)
FetchContent_MakeAvailable(or-tools)

# pybind11
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v3.0.1
)
set(PYBIND11_FINDPYTHON ON)
FetchContent_MakeAvailable(pybind11)

# Core library with implementation
add_library(sirius_core STATIC
  src/utils.cpp
  src/tables.cpp
  src/config.cpp
  src/solver.cpp
  src/logging.cpp
  src/instance.cpp
  src/interface.cpp
  src/cli.cpp
)
target_include_directories(sirius_core PUBLIC "${CMAKE_CURRENT_LIST_DIR}/include")
target_link_libraries(sirius_core PRIVATE ortools::ortools Threads::Threads)

# Native executable for the wrapper to run
add_executable(sirius src/sirius_main.cpp)   # tiny main that calls sirius_main
target_link_libraries(sirius PRIVATE sirius_core)

# Python extension: _sirius (produces _sirius.cpython-*.so with proper PyInit)
pybind11_add_module(_sirius src/pybind_module.cpp)

# target_link_libraries(_sirius PRIVATE sirius_core)

# make sure the .so lands in ${CMAKE_BINARY_DIR}
set_target_properties(_sirius PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
