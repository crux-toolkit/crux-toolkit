project(crux)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

set (CRUX_VERSION_MAJOR "1")
set (CRUX_VERSION_MINOR "40")
set (CRUX_VERSION ${CRUX_VERSION_MAJOR}.${CRUX_VERSION_MINOR})

add_subdirectory(doc/user)
add_subdirectory(src/external)
add_subdirectory(src/c)
add_subdirectory(src/c/test/smoke-tests)
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

# Set default TARGET_ARCH
if (NOT TARGET_ARCH)
 set(TARGET_ARCH "${CMAKE_SYSTEM}")
endif (NOT TARGET_ARCH)

set (CPACK_GENERATOR ZIP)
set (CPACK_SOURCE_GENERATOR TGZ)
set (CPACK_PACKAGE_VERSION_MAJOR ${CRUX_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${CRUX_VERSION_MINOR})
set (CRUX_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR})
set(CPACK_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CMAKE_SYSTEM_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.source")
set (
 CPACK_SOURCE_IGNORE_FILES 
   /_CPack_Packages/;
   /CMakeCache.txt/;
   /CMakeFiles/;
   /\\\\.svn/;
   /\\\\.svnignore/;
   /\\\\.swp/;
   /abstract/;
   /asms2007/;
   /benchmark/;
   /bin/;
   /build/;
   /crux.sln/;
   /developer/;
   /download/;
   /guispec/;
   /mods/;
   /open_source_download/;
   /paper/;
   /percolator/;
   /results/;
   /test/;
   /windows/;
   /*\\\\.cmake$/;
   /*\\\\.gz$/
   /*\\\\.zip$/
)
include (CPack)
