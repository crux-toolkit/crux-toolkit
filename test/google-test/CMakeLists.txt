 # If GoogleTest is not installed, build and install a private copy.
# Skip this step if GOOGLE_TEST is not defined.
if(GOOGLE_TEST)
    message(STATUS "Downloading GoogleTest package.")
    configure_file(GTestExtProject.cmake.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
      message(FATAL_ERROR "Failed to download GoogleTest: ${result}.")
    endif()
    message(STATUS "Building GoogleTest.")
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
      message(FATAL_ERROR "Failed to build GoogleTest package: ${result}.")
    endif()
    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                     EXCLUDE_FROM_ALL)
    set(GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include")
    message(STATUS "GoogleTest has been installed succesfully: ${GTEST_INCLUDE_DIRS} ${GTEST_LIBRARY_PATH}.")
endif(GOOGLE_TEST)

# Linking and building unit tests
include_directories(${GTEST_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_SOURCE_DIR}/src)

link_directories(
    ${CMAKE_BINARY_DIR}/app/bullseye
    ${CMAKE_BINARY_DIR}/app/hardklor
    ${CMAKE_BINARY_DIR}/ext/lib
    ${CMAKE_BINARY_DIR}/ext/lib64
    ${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src
    ${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/blas
    ${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/fido
    ${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/picked_protein)
    
add_executable(gtest_unit
    UnitTest_crux_main.cpp
    UnitTest_GoogleAnalytics.cpp)

# Flags for generating coverage data
# if(COVERAGE)
#   target_compile_options(perclibrary PUBLIC -ftest-coverage -fprofile-arcs)
#   target_compile_options(gtest_unit PUBLIC -ftest-coverage -fprofile-arcs)
#   target_link_libraries(gtest_unit -fprofile-arcs)
# endif(COVERAGE)
target_link_libraries(gtest_unit gtest gtest_main
    bullseye
    crux-support
    tide-support
    cometsearch
    kojaksearch
    gflags
    hardklor
    mstoolkitlite
    perclibrary
    blas
    fido
    picked_protein
    protobuf
    dl
    m
    pthread
    pwiz_data_identdata
    pwiz_data_identdata_version
    pwiz_data_misc
    pwiz_data_msdata
    pwiz_data_msdata_mz5
    pwiz_data_msdata_version
    pwiz_data_proteome
    pwiz_data_common
    pwiz_utility_chemistry
    pwiz_utility_minimxml
    pwiz_utility_misc
    pwiz_version
    boost_chrono
    boost_filesystem
    boost_iostreams
    boost_serialization
    boost_system
    boost_thread
    hdf5pp
    hdf5
    ${RT_LIBRARY}
    SHA1
    zlib
  )
  add_test(UnitTest_crux_RunAllTests gtest_unit)

# Important to use relative paths here (used by CPack)!
#install(TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin)