cmake_minimum_required(VERSION 3.6)
project(kojak)
# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(
  ${BASE_INCLUDE_PATH} 
  ${MSTOOLKIT_INCLUDE_PATH} 
  ${HARDKLOR_INCLUDE_PATH}
  ${CMAKE_CURRENT_BINARY_DIR}
)

#############################################################################
# COMPILE Kojak
###############################################################################

set(
  KOJAK_SOURCE
  KojakManager.cpp 
  KParams.cpp 
  KAnalysis.cpp 
  KData.cpp 
  KDB.cpp 
  KPrecursor.cpp 
  KSpectrum.cpp 
  KIons.cpp 
  KIonSet.cpp 
  KLog.cpp 
  KTopPeps.cpp 
  Threading.cpp 
  CometDecoys.cpp
)
add_library(
  kojaksearch 
  STATIC 
  ${KOJAK_SOURCE}
)

if (UNIX)
  target_compile_options(
    kojaksearch
    PUBLIC
    -O3 -static -std=c++11 -I. -I${HEADER_PATH} -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC
  )
elseif(WIN32 AND NOT CYGWIN)
  target_compile_definitions(
    kojaksearch
    PUBLIC
    _CRT_SECURE_NO_WARNINGS
  )
endif(UNIX)
  

################################################################################
								  
# INSTALL Kojak
install(TARGETS kojaksearch EXPORT kojaksearch DESTINATION lib)
install(DIRECTORY . DESTINATION include/kojak FILES_MATCHING PATTERN "*.h") 


