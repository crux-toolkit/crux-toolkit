# toolchain file for building a 32bit version on a 64bit host

# use it like this:
# cmake -DCMAKE_TOOLCHAIN_FILE=linux_i686.toolchain.cmake <sourcedir>

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_VERSION 1)
set(CMAKE_SYSTEM_PROCESSOR "i686")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32" CACHE STRING "c++ flags")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -m32" CACHE STRING "c flags")
set(CMAKE_LD_FLAGS   "${CMAKE_LD_FLAGS} -m32" CACHE STRING "linker flags")
set(BUILD_32 ON)
