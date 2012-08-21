include(ExternalProject)
ExternalProject_Add(
  XercesC
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/xerces
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/xerces/configure --disable-net --prefix=${CMAKE_CURRENT_SOURCE_DIR}/xerces
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/xerces
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)
ExternalProject_Add(
  percolator
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/percolator
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/percolator
)
add_dependencies(percolator XercesC)
add_subdirectory(MSToolkit)

