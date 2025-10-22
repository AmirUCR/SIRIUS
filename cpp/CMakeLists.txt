cmake_minimum_required(VERSION 3.20)
project(SIRIUS LANGUAGES CXX)

if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  add_compile_options(-march=x86-64 -mtune=generic)
endif()

# --- Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Threads (robust)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads QUIET)
if(Threads_FOUND)
  add_library(pthread_if INTERFACE)
  target_link_libraries(pthread_if INTERFACE Threads::Threads)
else()
  add_library(pthread_if INTERFACE)
  target_compile_options(pthread_if INTERFACE -pthread)
  target_link_options(pthread_if INTERFACE -pthread)
endif()

# --- OR-Tools: vendored
set(OR_TOOLS_DIR "${CMAKE_CURRENT_LIST_DIR}/ortools")
if(NOT EXISTS "${OR_TOOLS_DIR}")
  message(FATAL_ERROR "Vendored OR-Tools directory not found: ${OR_TOOLS_DIR}")
endif()

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

# --- Includes
target_include_directories(sirius PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/include"
  "${OR_TOOLS_DIR}/include"
)

# --- RPATH so the wheel’s sirius/lib is searched at runtime
set(_rpath "\$ORIGIN/../lib")
set_target_properties(sirius PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  BUILD_RPATH   "${_rpath}"
  INSTALL_RPATH "${_rpath}"
)

# --- Link order consolidated (m, dl, pthreads + OR-Tools)
#   Use 'm' and 'dl' (CMake adds -l), and the pthread_if interface
target_link_libraries(sirius PRIVATE
  "${ORTOOLS_LIB}"
  m
  dl
  pthread_if
)

# --- Release opts
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sirius PRIVATE -O3 -pipe)
  endif()
endif()
