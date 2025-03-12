cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW)
project(cometsearch)
# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(
  ${BASE_INCLUDE_PATH} 
  ${MSTOOLKIT_INCLUDE_PATH} 
  ${CMAKE_CURRENT_BINARY_DIR}/CometSearch
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

#############################################################################
# COMPILE MSToolkit
###############################################################################

set(
  COMETSEARCH
  CometSearch/Threading.cpp 
  CometSearch/CombinatoricsUtils.cpp
  CometSearch/CometFragmentIndex.cpp
  CometSearch/CometInterfaces.cpp 
  CometSearch/CometSearch.cpp 
  CometSearch/CometPeptideIndex.cpp 
  CometSearch/CometPreprocess.cpp 
  CometSearch/CometPostAnalysis.cpp 
  CometSearch/CometMassSpecUtils.cpp 
  CometSearch/CometWriteMzIdentML.cpp
  CometSearch/CometWriteOut.cpp 
  CometSearch/CometWriteSqt.cpp 
  CometSearch/CometWritePepXML.cpp 
  CometSearch/CometWritePercolator.cpp 
  CometSearch/CometWriteTxt.cpp 
  CometSearch/CometSearchManager.cpp 
  CometSearch/ModificationsPermuter.cpp 
)
 
add_library(
  cometsearch 
  STATIC 
  ${COMETSEARCH}
)
								  
if (UNIX)
  if (CMAKE_BUILD_TYPE MATCHES "Release")
    target_compile_options  (
      cometsearch
      PUBLIC
      -O3 -fpermissive -Wall -Wextra -Wno-write-strings -static -std=c++14 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DCRUX -I. -I$(MSTOOLKIT)/include
    )
  else()
    target_compile_options  (
      cometsearch
      PUBLIC
      -g -O0 -fpermissive -Wall -Wextra -Wno-write-strings -static -std=c++14 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC  -DCRUX -I. -I$(MSTOOLKIT)/include
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
  

################################################################################
								  
# INSTALL cometsearch library
install(TARGETS cometsearch EXPORT cometsearch DESTINATION lib)
install(DIRECTORY CometSearch/ DESTINATION include/CometSearch FILES_MATCHING PATTERN "*.h") 


