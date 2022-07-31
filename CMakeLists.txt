cmake_minimum_required(VERSION 3.0.0)
project(renderer VERSION 0.1.0)

include(CTest)
enable_testing()
SET(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_CXX_STANDARD 17)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
include_directories(include)

add_executable(renderer src/main.cpp src/model.cpp src/tgaimage.cpp)
