link_directories(${PERCOLATOR_BINARY_DIR}/src)

file(GLOB QVALITY_SOURCES *.cpp)

add_executable(qvality ${QVALITY_SOURCES})

target_link_libraries(qvality perclibrary fido)

# Relative paths don't seem to work on Apple
install(TARGETS qvality DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
