project(crux)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

add_subdirectory(doc/user)
add_subdirectory(src/external/MSToolkit)
add_subdirectory(src/external/)
add_subdirectory(src/c)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)

set (CPACK_SOURCE_GENERATOR TGZ)
set (CPACK_PACKAGE_VERSION_MAJOR 1)
set (CPACK_PACKAGE_VERSION_MINOR 8)
set (CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES};/abstract/;/asms2007/;/developer/;/guispec/;/mods/;/paper/)
include (CPack)
