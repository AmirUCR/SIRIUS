cmake_minimum_required(VERSION 3.20)
project(SIRIUS LANGUAGES CXX)

# --- Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- OR-Tools: local/vendored installation
set(OR_TOOLS_DIR "${CMAKE_CURRENT_LIST_DIR}/ortools")

list(APPEND CMAKE_PREFIX_PATH "${OR_TOOLS_DIR}/lib/cmake" "${OR_TOOLS_DIR}/share/eigen3/cmake/")
find_package(ZLIB NO_MODULE REQUIRED)
find_package(ortools NO_MODULE REQUIRED)

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

target_include_directories(sirius PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/include"
  "${OR_TOOLS_DIR}/include"
)

# --- Linking
target_link_libraries(sirius PRIVATE
  ortools::ortools
  -pthread
)

# Emit binary into the build dir (setup.py copies build_temp/sirius)
set_target_properties(sirius PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)