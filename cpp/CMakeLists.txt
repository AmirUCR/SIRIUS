# cpp/CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(SIRIUS LANGUAGES CXX)

# --- Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

# --- OR-Tools: vendored
set(OR_TOOLS_DIR "${CMAKE_CURRENT_LIST_DIR}/ortools")

if(NOT EXISTS "${OR_TOOLS_DIR}")
  message(FATAL_ERROR "Vendored OR-Tools directory not found: ${OR_TOOLS_DIR}")
endif()

# Find the library file
find_library(ORTOOLS_LIB
  NAMES ortools
  PATHS "${OR_TOOLS_DIR}/lib"
  NO_DEFAULT_PATH
)

if(NOT ORTOOLS_LIB)
  message(FATAL_ERROR "libortools (.so or .a) not found in ${OR_TOOLS_DIR}/lib")
endif()

# --- Sources
add_executable(sirius
  sirius.cpp
  src/utils.cpp
  src/tables.cpp
  src/config.cpp
  src/solver.cpp
  src/logging.cpp
  src/instance.cpp
  src/interface.cpp
)

# Headers live one level up: <repo_root>/sirius/include
target_include_directories(sirius PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/include"
  "${OR_TOOLS_DIR}/include"
)

target_link_libraries(sirius PRIVATE
  "${ORTOOLS_LIB}"
  Threads::Threads
)

# Emit binary into the build dir (setup.py copies build_temp/sirius)
set_target_properties(sirius PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Optional: a touch of optimization for Release
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sirius PRIVATE -O3 -pipe)
  endif()
endif()
