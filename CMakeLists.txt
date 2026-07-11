# Alternative build via CMake (Visual Studio, Ninja, or GCC/Clang for the test)
#   mkdir build && cd build
#   cmake -A x64 ..          (Windows / Visual Studio)
#   cmake --build . --config Release
cmake_minimum_required(VERSION 3.15)
project(CorelNestEngine CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(CorelNestEngine SHARED CorelNestEngine.cpp)
if(WIN32)
  target_sources(CorelNestEngine PRIVATE CorelNestEngine.def)
  # static CRT -> single-file DLL with no runtime dependency
  set_property(TARGET CorelNestEngine PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded")
endif()

add_executable(nest_test test_harness.cpp CorelNestEngine.cpp)
enable_testing()
add_test(NAME engine_sanity COMMAND nest_test)
