cmake_minimum_required(VERSION 3.10)
project(SIRIUS)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)

# Set the path to the OR-Tools library relative to the project directory
set(OR_TOOLS_DIR "${CMAKE_SOURCE_DIR}/ortools")

# Add OR-Tools include directories and libraries
include_directories(${OR_TOOLS_DIR}/include)
link_directories(${OR_TOOLS_DIR}/lib)

# Add the executable
add_executable(sirius sirius.cpp)

# Link the OR-Tools library
target_link_libraries(sirius ortools)
