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
  WalkIntriangulation ${CMAKE_CURRENT_SOURCE_DIR}/Data/StandardTestMesh.vtk 1 2 0
  )

