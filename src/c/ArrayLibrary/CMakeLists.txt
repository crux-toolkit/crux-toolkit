project(ArrayLibrary)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

add_library(
  ArrayLibrary
  STATIC
  MatrixFunctions.cpp
  Numerical.cpp
  PackedMatrix.cpp
  PackedSetFunctions.cpp
  PackedVec.cpp
  Random.cpp
  VectorFunctions.cpp
)
