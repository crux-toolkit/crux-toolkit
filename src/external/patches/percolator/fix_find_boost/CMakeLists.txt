###############################################################################
# SET COMPILER OPTIONS
###############################################################################

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_C_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
IF(APPLE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  
  set(LINK_FLAGS "-mmacosx-version-min=10.5")
ENDIF(APPLE)

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

#find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
#if(XERCESC_FOUND)
#  message(STATUS "Package Xerces-C++ found: ${XERCESC_INCLUDE_DIR}")
#else(XERCESC_FOUND)
#  message(FATAL_ERROR "The package Apache Xerces-C++ has not been found")
#endif(XERCESC_FOUND)
#include_directories(${XERCESC_INCLUDE_DIR})

# Set Boost and Xerces variables
include_directories(${CMAKE_PREFIX_PATH}/include)
link_directories(${CMAKE_PREFIX_PATH}/lib)
if (WIN32 AND NOT CYGWIN)
  set(
    Boost_LIBRARIES 
    ${CMAKE_PREFIX_PATH}/lib/libboost_filesystem-vc100-mt.lib
    ${CMAKE_PREFIX_PATH}/lib/libboost_system.-vc100-mt.lib
  )
  set(
    XERCESC_LIBRARIES 
    ${CMAKE_PREFIX_PATH}/lib/xerces-c_3D.lib 
  )
else (WIN32 AND NOT CYGWIN)
  set(
    Boost_LIBRARIES 
    ${CMAKE_PREFIX_PATH}/lib/libboost_filesystem.a 
    ${CMAKE_PREFIX_PATH}/lib/libboost_system.a
  )
  set(
    XERCESC_LIBRARIES 
    ${CMAKE_PREFIX_PATH}/lib/libxerces-c.a
  )
endif (WIN32 AND NOT CYGWIN)

set(Boost_LIBRARIES ${CMAKE_PREFIX_PATH}/lib/libboost_filesystem.a ${CMAKE_PREFIX_PATH}/lib/libboost_system.a)
set(XERCES_LIBRARIES ${CMAKE_PREFIX_PATH}/lib/ ${CMAKE_PREFIX_PATH}/lib/libboost_system.a)

if (UNIX)
  find_package(Pthreads)
  if(PTHREADS_FOUND)
    message(STATUS  "Pthreads found")
  else(PTHREADS_FOUND)
    message(FATAL_ERROR "The package Pthreads has not been found")
  endif()
endif (UNIX)

###############################################################################
# RUN CODESYNTHESIS
###############################################################################

set(xsdfiles_in percolator_in.cxx percolator_in.hxx )
set(percolator_in xml/percolator_in.xsd)
if (CYGWIN)
    # On Cygwin XSD expects Windows style paths so we use cygpath to covert from UNIX style paths.
	message(STATUS "Running on cygwin")
	execute_process(
		COMMAND cygpath -w ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
		OUTPUT_VARIABLE XSD_IN_INPUT_FILE
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else(CYGWIN)
	message(STATUS "Not running on Cygwin")
	set(XSD_IN_INPUT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in})
endif(CYGWIN)
message(STATUS "XSD_IN_INPUT FILE = ${XSD_IN_INPUT_FILE}")
message(STATUS "xsdfiles_in = ${xsdfiles_in}")

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pout (percolator output)
set(xsdfiles_out percolator_out.cxx percolator_out.hxx)
set(percolator_out xml/percolator_out.xsd)
if (CYGWIN)
    # On Cygwin XSD expects Windows style paths so we use cygpath to covert from UNIX style paths.
	message(STATUS "Running on Cygwin")
	execute_process(
		COMMAND cygpath -w ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_out}
		OUTPUT_VARIABLE XSD_OUT_INPUT_FILE
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else(CYGWIN)
	message(STATUS "Not running on Cygwin")
	set(XSD_OUT_INPUT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_out})
endif(CYGWIN)
message(STATUS "XSD_OUT_INPUT_FILE = ${XSD_OUT_INPUT_FILE}")
message(STATUS "xsdfiles_out = ${xsdfiles_out}")

add_definitions(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)


###############################################################################
# COMPILE INTERNAL LIBRARIES
###############################################################################

#COMPILE FIDO
add_subdirectory(fido)

#INCLUDE PERCOLATOR HEADERS
include_directories(.) 

#INCLUDE FIDO HEADERS FOR PERCLIBRARY 
include_directories(. fido)


#############################################################################
# PERCOLATOR
###############################################################################

add_library(perclibrary STATIC ${xsdfiles_in} ${xsdfiles_out} parser.cxx serializer.cxx BaseSpline.cpp DescriptionOfCorrect.cpp MassHandler.cpp PSMDescription.cpp ResultHolder.cpp SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp Enzyme.cpp Globals.cpp Normalizer.cpp PercolatorCInterface.cpp SanityCheck.cpp UniNormalizer.cpp DataSet.cpp FeatureNames.cpp LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp ProteinProbEstimator.cpp Scores.cpp SqtSanityCheck.cpp ssl.cpp EludeModel.cpp)

# COMPILE PERCOLATOR
add_executable(percolator main.cpp)

IF(MINGW OR WIN32)
  set_property(TARGET percolator PROPERTY LINK_SEARCH_END_STATIC TRUE)
ENDIF()

if(MINGW OR WIN32)
  set_target_properties(percolator PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
endif()


 
# SET LIBRARIES FOR PERCOLATOR
if (UNIX)
  target_link_libraries(
    percolator 
    perclibrary 
    fido 
    pthread 
    ${PTHREADS_LIBRARIES} 
    ${XERCESC_LIBRARIES} 
    ${Boost_LIBRARIES}
  )
else (UNIX)
  target_link_libraries(
    percolator 
    perclibrary 
    fido 
    ${XERCESC_LIBRARIES} 
    ${Boost_LIBRARIES}
  )
endif (UNIX)
  
 

# INSTALL PERCOLATOR
IF(NOT MINGW AND NOT WIN32)
  install(TARGETS percolator EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
endif()

###############################################################################
# COMPILE QUALITY
###############################################################################

add_subdirectory(qvality)
