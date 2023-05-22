# CMakeList.txt : CMake project for crux-quant, include source and define
# project specific logic here.
#
cmake_minimum_required(VERSION 3.15)

project ("crux")

set (CRUX_QUANT_VERSION_MAJOR "0")
set (CRUX_QUANT_VERSION_MINOR "0")
set (CRUX_QUANT_VERSION ${CRUX_QUANT_VERSION_MAJOR}.${CRUX_QUANT_VERSION_MINOR})


macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  else (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_STANDARD 11)
  endif (CMAKE_VERSION VERSION_LESS "3.1")
endmacro(use_cxx11)

use_cxx11()


enable_testing()

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(crux_test crux_quant_test.cpp)

target_link_libraries(crux_test crux-quant-support)
target_link_libraries(crux_test  GTest::gtest GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(crux_test)