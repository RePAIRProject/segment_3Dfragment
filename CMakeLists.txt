cmake_minimum_required(VERSION 3.16)
project(segment_3Dfragment)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
option(LIBIGL_PNG "Build target igl::png" ON)
include(libigl)

# Add your project files
file(GLOB SRC_FILES src/*.cpp *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})

# Ensure the necessary include directories
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Link against required libraries
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl::png)
