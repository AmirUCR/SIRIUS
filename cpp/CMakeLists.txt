cmake_minimum_required(VERSION 3.24)
project(SIRIUS VERSION 1.0 LANGUAGES CXX)

# --- Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(FetchContent)
FetchContent_Declare(
  or-tools
  GIT_REPOSITORY https://github.com/google/or-tools.git
  GIT_TAG        "v9.14"
)

# After the following call, the CMake targets defined by OR-Tools
# will be defined and available to the rest of the build
FetchContent_MakeAvailable(or-tools)

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