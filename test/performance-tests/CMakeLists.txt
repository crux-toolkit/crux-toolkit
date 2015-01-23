cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

# Copy files for performance tests in case we are building out of source.
file(COPY crux.param DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY stored-performance.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY worm+contaminants.fa DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY performance.html DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY 051708-worm-ASMS-10.ms2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY run-performance-test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY runall DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY stored-xcorr.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(performance-tests COMMAND ./run-performance-test.sh)
