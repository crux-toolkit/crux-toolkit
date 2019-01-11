project(crux)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

if (BUILD_32)
  message("-- 32 bit build selected")
endif(BUILD_32)


set (CRUX_VERSION_MAJOR "3")
set (CRUX_VERSION_MINOR "2")
set (CRUX_VERSION ${CRUX_VERSION_MAJOR}.${CRUX_VERSION_MINOR})

message(
  "-- Creating build files for Crux ${CRUX_VERSION} for system ${CMAKE_SYSTEM} " 
  "(short name ${CMAKE_SYSTEM_NAME})"
)

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

add_subdirectory(doc/user)
add_subdirectory(ext)
add_subdirectory(src)


# Only process smoke tests directory if it exiists
if (EXISTS "${CMAKE_SOURCE_DIR}/test/smoke-tests/runall")
  add_subdirectory(test/smoke-tests)
endif (EXISTS "${CMAKE_SOURCE_DIR}/test/smoke-tests/runall")

# Only process performance tests directory if it exiists
if (EXISTS "${CMAKE_SOURCE_DIR}/test/performance-tests/runall")
  add_subdirectory(test/performance-tests)
endif (EXISTS "${CMAKE_SOURCE_DIR}/test/performance-tests/runall")

install(FILES license.txt DESTINATION .)

add_custom_target(
  deep-clean
  COMMAND ${CMAKE_COMMAND} --build . --target clean
  COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/build
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/downloads
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/CometSearch
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/MSToolkit
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/boost
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/gflags-2.1.2
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/google
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/pwiz
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/xercesc
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/include/xsd
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/lib/pkgconfig
  COMMAND ${CMAKE_COMMAND} -E remove ext/lib/*.a
  COMMAND ${CMAKE_COMMAND} -E remove ext/lib/*.la
  COMMAND ${CMAKE_COMMAND} -E remove ext/lib/*.lib
  COMMAND ${CMAKE_COMMAND} -E remove_directory ext/share
  COMMENT "Cleaning project tree."
  COMMENT "Deleting all files for external tools."
)

add_custom_target(
  externals-only 
  COMMAND ${CMAKE_MAKE_PROGRAM}
  WORKING_DIRECTORY ext
  COMMENT "Building only external dependendcies."
)

add_custom_target(
  no-externals 
  COMMAND ${CMAKE_MAKE_PROGRAM} crux  xlink-score-spectrum xlink-assign-ions 
  WORKING_DIRECTORY src
  COMMENT "Not building external dependendcies."
)

# a custom target for obtaining the SVN revision number
add_custom_target(crux_version ALL DEPENDS crux_version_header )

# creates crux_version.h using cmake script
add_custom_command(
    OUTPUT crux_version_header 
    ${CMAKE_CURRENT_BINARY_DIR}/crux_version.h
    COMMAND 
      ${CMAKE_COMMAND} 
        -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DCRUX_VERSION=${CRUX_VERSION}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/make-version-header.cmake
)

add_dependencies(barista comet)
add_dependencies(bullseye hardklor)
add_dependencies(comet ProteoWizard)
add_dependencies(crux barista)
add_dependencies(crux bullseye)
add_dependencies(crux crux-support)
add_dependencies(crux hardklor)
add_dependencies(crux percolator)
add_dependencies(crux tide-support)
add_dependencies(crux xlink)
add_dependencies(crux-support crux_version)
add_dependencies(crux-support percolator)
add_dependencies(hardklor comet)
add_dependencies(percolator comet)
add_dependencies(tide-support comet)
add_dependencies(tide-support protobuf-3.0.0)
add_dependencies(tide-support gflags-2.2.1)
add_dependencies(xlink comet)
add_dependencies(xlink percolator)
add_dependencies(xlink ProteoWizard)
add_dependencies(create-docs crux)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)

if (WIN32 AND NOT CYGWIN)
  # If windows build we need to select 32/64 bit platform
  if (CMAKE_GENERATOR MATCHES "^.*Win64$")
    SET(BUILD_PROCESSOR "AMD64")
  else (CMAKE_GENERATOR MATCHES "^.*Win64$")
    SET(BUILD_PROCESSOR "i386")
  endif (CMAKE_GENERATOR MATCHES "^.*Win64$")
else (WIN32 AND NOT CYGWIN)
  SET(BUILD_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR})
endif (WIN32 AND NOT CYGWIN)

set (CPACK_GENERATOR ZIP)
set (CPACK_SOURCE_GENERATOR TGZ)
set (CPACK_PACKAGE_VERSION_MAJOR ${CRUX_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${CRUX_VERSION_MINOR})
set (CRUX_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR})
set(CPACK_PACKAGE_FILE_NAME
  ${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CMAKE_SYSTEM_NAME}.${BUILD_PROCESSOR})
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.Source")
set (
 CPACK_SOURCE_IGNORE_FILES 
   /\\\\..*$;
   \\\\./_CPack_Packages/;
   /benchmark/;
   /bin/;
   /results/;
   /doc/developer/;
   /doc/guispec/;
   /doc/mods/;
   /doc/user/download/;
   /doc/user/open_source_download/;
   /test/;
   /src/CruxGui/;
   /ext/bin/;
   /ext/build/;
   /ext/downloads/;
   /ext/include/;
   /ext/lib/;
   /ext/share/;
   /src/python/;
   /CMakeCache.txt$;
   /CMakeFiles/;
   \\\\.gz$;
   \\\\.swp$;
   \\\\.zip$;
   ${CPACK_SOURCE_IGNORE_FILES}
)

# Generate checksum files for binary and source packages.
add_custom_target(
  checksum
 	COMMAND ${CMAKE_COMMAND} -E  md5sum ${CPACK_PACKAGE_FILE_NAME}.zip > ${CPACK_PACKAGE_FILE_NAME}.zip.md5sum.txt
 	COMMAND ${CMAKE_COMMAND} -E  md5sum ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz > ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz.md5sum.txt
  COMMENT "Creating checksum files."
)

include (CPack)
