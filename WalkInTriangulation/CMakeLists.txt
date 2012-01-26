cmake_minimum_required(VERSION 2.8.4)

project( WalkInTriangulation )

FIND_PACKAGE ( ITK REQUIRED) 
INCLUDE( ${USE_ITK_FILE} ) 

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src )

# create the executable that will be used by all tests
ADD_EXECUTABLE( WalkInTriangulation
  src/itkWalkInTriangulationFunctor.h
  src/WalkInTriangulation.cxx
  src/predicates.c
  )
TARGET_LINK_LIBRARIES( WalkInTriangulation
${ITK_LIBRARIES}
  )

INCLUDE( CTest )

# set of test in order to validate the algorithm
# validation is done on a 5x5 points square tile triangulated mesh
# we check the path taken in 6 differents cases
# 
ADD_TEST( WalkInTriangulationTest01
  WalkIntriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 3.75 3.25 0 7 0 3 10 13 20 23 30
  )
ADD_TEST( WalkInTriangulationTest02
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.25 0.75 31 7 31 28 21 18 11 8 1
  )
ADD_TEST( WalkInTriangulationTest03
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.25 3.75 0 8 0 1 8 9 16 17 24 25
  )
ADD_TEST( WalkInTriangulationTest04
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.75 0.25 25 8 25 24 17 16 9 8 1 0
  )
ADD_TEST( WalkInTriangulationTest05
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 3.75 0.25 0 7 0 3 2 5 4 7 6
  )
ADD_TEST( WalkInTriangulationTest06
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.75 0.25 6 7 6 7 4 5 2 3 0
  )
