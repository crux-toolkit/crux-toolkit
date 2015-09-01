cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)
install(
  DIRECTORY . 
  DESTINATION doc/comet-options
  PATTERN "CMakeLists.txt" EXCLUDE
)
