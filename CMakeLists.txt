project(crux)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

set (CRUX_VERSION_MAJOR "1")
set (CRUX_VERSION_MINOR "40")
set (CRUX_VERSION ${CRUX_VERSION_MAJOR}.${CRUX_VERSION_MINOR})

add_subdirectory(doc/user)
add_subdirectory(src/external)
add_subdirectory(src/c)

# Only process smoke tests directory if it exiists
if (EXISTS "${CMAKE_SOURCE_DIR}/src/c/test/smoke-tests/runall")
  add_subdirectory(src/c/test/smoke-tests)
endif (EXISTS "${CMAKE_SOURCE_DIR}/src/c/test/smoke-tests/runall")

install(FILES license.txt DESTINATION .)

add_custom_target(
  deep-clean
  COMMAND ${CMAKE_MAKE_PROGRAM} clean
  COMMAND ${CMAKE_MAKE_PROGRAM} deep-clean-external
  COMMAND find ${CMAKE_BINARY_DIR} -name CMakeCache.txt -exec rm "{}" "\;"
  COMMENT "Cleaning project tree."
)

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
   \\\\.cmake$;
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
