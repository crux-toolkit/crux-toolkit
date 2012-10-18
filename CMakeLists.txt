project(crux)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

add_subdirectory(doc/user)
add_subdirectory(src/external/MSToolkit)
add_subdirectory(src/external/)
add_subdirectory(src/c)

install (FILES license.txt DESTINATION .)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)

set (CPACK_GENERATOR ZIP)
set (CPACK_SOURCE_GENERATOR TGZ)
set (CPACK_PACKAGE_VERSION_MAJOR "1")
set (CPACK_PACKAGE_VERSION_MINOR "39")
set (CRUX_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR})

configure_file (
 "src/c/CruxConfig.h.in"
 "src/c/CruxConfig.h"
 )

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
   /results/;
   /test/;
   /windows/;
   /*\\\\.cmake$/;
   /*\\\\.gz$/
)
include (CPack)
