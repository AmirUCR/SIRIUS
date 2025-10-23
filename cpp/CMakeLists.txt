cmake_minimum_required(VERSION 3.20)
project(SIRIUS LANGUAGES CXX)

# --- Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- OR-Tools: local/vendored installation
set(OR_TOOLS_DIR "${CMAKE_CURRENT_LIST_DIR}/ortools")

if(NOT EXISTS "${OR_TOOLS_DIR}")
  message(FATAL_ERROR "Vendored OR-Tools directory not found: ${OR_TOOLS_DIR}")
endif()

# Guide CMake to find the OR-Tools configuration files within the vendored directory.
# This replaces the manual 'find_library' and enables correct transitive dependency linking.
list(APPEND CMAKE_PREFIX_PATH "${OR_TOOLS_DIR}")
find_package(OR_TOOLS REQUIRED)

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
  "${OR_TOOLS_DIR}/include" # Keep this to ensure OR-Tools headers are found
)

# --- Linking
target_link_libraries(sirius PRIVATE
  # Link to the main OR-Tools CMake target. This target is guaranteed to include
  # all required transitive libraries (like absl and the SAT solver's specific
  # definition files) resolving the previous 'undefined reference' errors.
  ortools::ortools
  
  # The -pthread flag is necessary for thread-safe applications (common with OR-Tools/Abseil)
  # and is kept here as requested, though CMake's 'Threads::Threads' is often preferred.
  -pthread
)

# Emit binary into the build dir (setup.py copies build_temp/sirius)
set_target_properties(sirius PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)