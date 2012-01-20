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

include( CTest )

add_test( WalkInTriangulationTest01
  WalkIntriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 3.75 3.25 0 7 0 3 10 13 20 23 30
  )
add_test( WalkInTriangulationTest02
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.25 0.75 31 7 31 28 21 18 11 8 1
  )
add_test( WalkInTriangulationTest03
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.25 3.75 0 8 0 1 8 9 16 17 24 25
  )
add_test( WalkInTriangulationTest04
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.75 0.25 25 8 25 24 17 16 9 8 1 0
  )
add_test( WalkInTriangulationTest05
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 3.75 0.25 0 7 0 3 2 5 4 7 6
  )
add_test( WalkInTriangulationTest06
  WalkInTriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 0.75 0.25 6 7 6 7 4 5 2 3 0
  )
