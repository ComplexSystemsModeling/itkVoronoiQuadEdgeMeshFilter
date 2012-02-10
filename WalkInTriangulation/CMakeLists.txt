cmake_minimum_required( VERSION 2.8 )

project( WalkInTriangulation )

find_package( ITK REQUIRED ) 
include( ${USE_ITK_FILE} ) 

include_directories( ${PROJECT_SOURCE_DIR}/src )
include_directories( ${PointInCircle_DIR} )

# add source file (.c .cc .cpp .cxx ...)
#
set( WalkInTriangulation_src
  src/WalkInTriangulation.cxx
  ${predicates.c_location}
  )

# add header file (.h .hxx .txx ...)
#
set( WalkInTriangulation_header
  src/itkWalkInTriangulationFunction.h
  )

# if generat by Xcode add headers to executable
#
if( ${CMAKE_GENERATOR} MATCHES "Xcode" )
  set( WalkInTriangulation_src 
    ${WalkInTriangulation_header}
    ${WalkInTriangulation_src}
    )
endif()

# create the executable that will be used by all tests
#
add_executable( WalkInTriangulation
  ${WalkInTriangulation_src}
  )

# add libraries
#
target_link_libraries( WalkInTriangulation
  ${ITK_LIBRARIES}
  )

# set of test in order to validate the algorithm
# validation is done on a 5x5 points square tile triangulated mesh
# we check the path taken in 6 differents cases
# 
add_test( WalkInTriangulationTest01
  ${CMAKE_CURRENT_BINARY_DIR}/WalkIntriangulation 3.75 3.25 0 7 0 3 10 13 20 23 30
  )
add_test( WalkInTriangulationTest02
  ${CMAKE_CURRENT_BINARY_DIR}/WalkInTriangulation 0.25 0.75 31 7 31 28 21 18 11 8 1
  )
add_test( WalkInTriangulationTest03
  ${CMAKE_CURRENT_BINARY_DIR}/WalkInTriangulation 0.25 3.75 0 8 0 1 8 9 16 17 24 25
  )
add_test( WalkInTriangulationTest04
  ${CMAKE_CURRENT_BINARY_DIR}/WalkInTriangulation 0.75 0.25 25 8 25 24 17 16 9 8 1 0
  )
add_test( WalkInTriangulationTest05
  ${CMAKE_CURRENT_BINARY_DIR}/WalkInTriangulation 3.75 0.25 0  7 0 3 2 5 4 7 6
  )
add_test( WalkInTriangulationTest06
  ${CMAKE_CURRENT_BINARY_DIR}/WalkInTriangulation 0.75 0.25 6 7 6 7 4 5 2 3 0
  )
