cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW)
project(cometsearch)

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(
  ${BASE_INCLUDE_PATH}
  ${MSTOOLKIT_INCLUDE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/CometSearch
  ${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/include
  ${CMAKE_CURRENT_BINARY_DIR}/CometSearch
  ${CMAKE_CURRENT_BINARY_DIR}/AScorePro/include
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

#############################################################################
# Build AScorePro static library (headers in AScorePro/include, sources in AScorePro/)
###############################################################################

# collect AScorePro sources; replace GLOB with explicit list if preferred
file(GLOB ASCOREPRO_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/*.c
  ${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/*.cc
)

if (ASCOREPRO_SRCS)
  add_library(
    ascorepro
    STATIC
    ${ASCOREPRO_SRCS}
  )

  # Force output base name to 'ascorepro' so on Unix we get libascorepro.a
#  set_target_properties(ascorepro PROPERTIES OUTPUT_NAME ascorepro)

  # Make headers available to dependents while building
  target_include_directories(
    ascorepro
    PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/include
      ${CMAKE_CURRENT_BINARY_DIR}/AScorePro/include
  )

  if (UNIX)
    if (CMAKE_BUILD_TYPE MATCHES "Release")
      target_compile_options(
        ascorepro
        PUBLIC
        -O3 -fpermissive -Wall -Wextra -Wno-write-strings -static -std=c++14 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DCRUX -I.
      )
    else()
      target_compile_options(
        ascorepro
        PUBLIC
        -g -O0 -fpermissive -Wall -Wextra -Wno-write-strings -static -std=c++14 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DCRUX -I.
      )
    endif()
  elseif(WIN32 AND NOT CYGWIN)
    target_compile_definitions(
      ascorepro
      PUBLIC
      _CRT_SECURE_NO_WARNINGS
      CRUX
    )
  endif()

  # Install AScorePro artifacts (headers and library) into the usual locations
  # ARCHIVE -> static libraries (.a) on Unix
  install(TARGETS ascorepro
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
  )

  # Copy the contents of AScorePro/include into include/AScorePro (no extra 'include' level)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/include/ DESTINATION include/AScorePro FILES_MATCHING PATTERN "*.h")
endif()

#############################################################################
# COMPILE COMETSEARCH
###############################################################################

set(
  COMETSEARCH
  CometSearch/CombinatoricsUtils.cpp
  CometSearch/CometAlignment.cpp
  CometSearch/CometFragmentIndex.cpp
  CometSearch/CometInterfaces.cpp
  CometSearch/CometMassSpecUtils.cpp
  CometSearch/CometModificationsPermuter.cpp
  CometSearch/CometPeptideIndex.cpp
  CometSearch/CometPostAnalysis.cpp
  CometSearch/CometPreprocess.cpp
  CometSearch/CometSearch.cpp
  CometSearch/CometSearchManager.cpp
  CometSearch/CometSpecLib.cpp
  CometSearch/CometWriteMzIdentML.cpp
  CometSearch/CometWritePepXML.cpp
  CometSearch/CometWritePercolator.cpp
  CometSearch/CometWriteSqt.cpp
  CometSearch/CometWriteTxt.cpp
  CometSearch/Threading.cpp
)

add_library(
  cometsearch
  STATIC
  ${COMETSEARCH}
)

# Make sure cometsearch can find AScorePro include files and MSToolkit headers
target_include_directories(
  cometsearch
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/CometSearch
    ${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/include
    ${CMAKE_CURRENT_BINARY_DIR}/AScorePro/include
    ${MSTOOLKIT_INCLUDE_PATH}
)

if (UNIX)
  if (CMAKE_BUILD_TYPE MATCHES "Release")
    target_compile_options(
      cometsearch
      PUBLIC
      -O3 -static -std=c++14 -fpermissive -Wall -Wextra -Wno-write-strings -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -D_NOSQLITE -D__int64=off64_t -I. -I${MSTOOLKIT_INCLUDE_PATH} -I${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/include
    )
  else()
    target_compile_options(
      cometsearch
      PUBLIC
      -O0 -static -std=c++14 -g -fpermissive -Wall -Wextra -Wno-write-strings -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -D_NOSQLITE -D__int64=off64_t -I. -I${MSTOOLKIT_INCLUDE_PATH} -I${CMAKE_CURRENT_SOURCE_DIR}/AScorePro/include
    )
  endif ()
elseif(WIN32 AND NOT CYGWIN)
  target_compile_definitions(
    cometsearch
    PUBLIC
    _CRT_SECURE_NO_WARNINGS
    CRUX
  )
endif(UNIX)

# Link cometsearch against AScorePro if that lib was created
if (TARGET ascorepro)
  target_link_libraries(cometsearch PUBLIC ascorepro)
endif()

################################################################################

# INSTALL cometsearch library and headers
install(TARGETS cometsearch
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CometSearch/ DESTINATION include/CometSearch FILES_MATCHING PATTERN "*.h")

# Note: AScorePro headers and static library are installed above with ascorepro target.
# With the OUTPUT_NAME set to 'ascorepro', the installed static library on Unix will be libascorepro.a
# and it will be placed in the same lib/ directory as libcometsearch.a
