cmake_minimum_required(VERSION 3.20)
project(SIRIUS LANGUAGES CXX)

# --- Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- OR-Tools: local/vendored installation
#set(OR_TOOLS_DIR "${CMAKE_CURRENT_LIST_DIR}/ortools")

#list(APPEND CMAKE_PREFIX_PATH "${OR_TOOLS_DIR}/lib/cmake" "${OR_TOOLS_DIR}/share/eigen3/cmake/")
#find_package(ZLIB NO_MODULE REQUIRED)
#find_package(ortools NO_MODULE REQUIRED)
#find_package(absl NO_MODULE REQUIRED)

include(FetchContent)
FetchContent_Declare(
  or-tools
  GIT_REPOSITORY https://github.com/google/or-tools.git
  GIT_TAG        "v9.11"
  # Passing multiple CMAKE variables to configure the OR-Tools build
  CMAKE_ARGS -DBUILD_DEPS:BOOL=ON
             -DBUILD_absl:BOOL=ON
             -DUSE_COINOR:BOOL=OFF
             -DUSE_HIGHS:BOOL=OFF
             -DUSE_SCIP:BOOL=OFF
             -DBUILD_SAMPLES:BOOL=OFF
             -DBUILD_EXAMPLES:BOOL=OFF
             -DUSE_DOTNET_8:BOOL=OFF
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

#target_include_directories(sirius PRIVATE
#  "${CMAKE_CURRENT_LIST_DIR}/include"
#  "${OR_TOOLS_DIR}/include"
#)

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