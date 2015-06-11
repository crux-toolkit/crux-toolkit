project(crux)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

set (CRUX_VERSION_MAJOR "2")
set (CRUX_VERSION_MINOR "1")
set (CRUX_VERSION ${CRUX_VERSION_MAJOR}.${CRUX_VERSION_MINOR})

message(
  "-- Creating build files for Crux ${CRUX_VERSION} for system ${CMAKE_SYSTEM} " 
  "(short name ${CMAKE_SYSTEM_NAME})"
)

add_subdirectory(doc/user)
add_subdirectory(src/external)
add_subdirectory(src/c)

# Only process smoke tests directory if it exiists
if (EXISTS "${CMAKE_SOURCE_DIR}/src/c/test/smoke-tests/runall")
  add_subdirectory(src/c/test/smoke-tests)
endif (EXISTS "${CMAKE_SOURCE_DIR}/src/c/test/smoke-tests/runall")

# Only process performance tests directory if it exiists
if (EXISTS "${CMAKE_SOURCE_DIR}/src/c/test/performance-tests/runall")
  add_subdirectory(src/c/test/performance-tests)
endif (EXISTS "${CMAKE_SOURCE_DIR}/src/c/test/performance-tests/runall")

install(FILES license.txt DESTINATION .)

add_custom_target(
  deep-clean
  COMMAND ${CMAKE_COMMAND} --build . --target clean
  COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/build
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/downloads
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/CometSearch
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/MSToolkit
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/boost
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/gflags
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/google
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/pwiz
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/xercesc
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/include/xsd
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/lib/pkgconfig
  COMMAND ${CMAKE_COMMAND} -E remove src/external/lib/*.a
  COMMAND ${CMAKE_COMMAND} -E remove src/external/lib/*.la
  COMMAND ${CMAKE_COMMAND} -E remove src/external/lib/*.lib
  COMMAND ${CMAKE_COMMAND} -E remove_directory src/external/share
  COMMENT "Cleaning project tree."
  COMMENT "Deleting all files for external tools."
)

add_custom_target(
  externals-only 
  COMMAND ${CMAKE_MAKE_PROGRAM}
  WORKING_DIRECTORY src/external
  COMMENT "Building only external dependendcies."
)

add_custom_target(
  no-externals 
  COMMAND ${CMAKE_MAKE_PROGRAM} crux  xlink-score-spectrum xlink-assign-ions 
  WORKING_DIRECTORY src/c
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
add_dependencies(crux posterior-error)
add_dependencies(crux tide-support)
add_dependencies(crux xlink)
add_dependencies(crux-support crux_version)
add_dependencies(crux-support percolator)
add_dependencies(hardklor comet)
add_dependencies(percolator ProteoWizard)
add_dependencies(tide-index tide-support)
add_dependencies(tide-search tide-support)
add_dependencies(tide-support ProtocolBuffers)
add_dependencies(tide-support google-gflags)
add_dependencies(xlink comet)
add_dependencies(xlink percolator)
add_dependencies(xlink ProteoWizard)
add_dependencies(xlink-assign-ions crux)
add_dependencies(xlink-score-spectrum crux)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)

execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
set (CPACK_GENERATOR ZIP)
set (CPACK_SOURCE_GENERATOR TGZ)
set (CPACK_PACKAGE_VERSION_MAJOR ${CRUX_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${CRUX_VERSION_MINOR})
set (CRUX_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR})
set(CPACK_PACKAGE_FILE_NAME
  ${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CMAKE_SYSTEM_NAME}.${ARCH})
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.Source")
set (
 CPACK_SOURCE_IGNORE_FILES 
   /\\\\..*$;
   \\\\./_CPack_Packages/;
   /benchmark/;
   /bin/;
   /results/;
   [^c]/tide/;
   /doc/developer/;
   /doc/guispec/;
   /doc/mods/;
   /doc/user/download/;
   /doc/user/open_source_download/;
   /src/c/test/;
   /src/CruxGui/;
   /src/external/bin/;
   /src/external/build/;
   /src/external/downloads/;
   /src/external/include/;
   /src/external/lib/;
   /src/external/share/;
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
add_custom_command(TARGET package POST_BUILD COMMAND echo "post package" COMMENT "Done with post package")
