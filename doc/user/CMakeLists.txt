cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

# Generate documentation
if (UNIX)
  set(CRUX_PATH "${PROJECT_BINARY_DIR}/src/crux")
elseif(WIN32)
  set(CMAKE_MSVCIDE_RUN_PATH ${PROJECT_BINARY_DIR}/ext/lib)
  set(CRUX_PATH "$<TARGET_FILE:crux>")
endif (UNIX)
add_custom_target(
  create-docs
  ALL
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc/user
  COMMAND 
    ${CMAKE_COMMAND} 
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
      -DCRUX_PATH:PATH=${CRUX_PATH}
      -DDOC_DIR:PATH=${PROJECT_SOURCE_DIR}/doc/user/commands
      -P  ${CMAKE_CURRENT_SOURCE_DIR}/create-docs.cmake
)

install(
  DIRECTORY . 
  DESTINATION doc
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "download" EXCLUDE 
)
