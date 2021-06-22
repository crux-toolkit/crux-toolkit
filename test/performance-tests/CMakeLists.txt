cmake_minimum_required(VERSION 3.6)
cmake_policy(VERSION 3.6)

# Copy files for performance tests in case we are building out of source.
file(COPY crux.param DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY worm+contaminants.fa DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY performance.html DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY 051708-worm-ASMS-10.ms2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY runall DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY run-performance-test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY runall DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY stored-plots DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(performance-tests COMMAND ./runall)
