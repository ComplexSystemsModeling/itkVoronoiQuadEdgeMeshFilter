cmake_minimum_required( VERSION 2.8 )

project( DelaunayTriangulation )

find_package( ITK REQUIRED) 
include( ${USE_ITK_FILE} ) 

include_directories( ${CMAKE_SOURCE_DIR}/src )
include_directories( ${PointInCircle_DIR} )
include_directories( ${WalkInTriangulation_DIR} )

# add source file
#
set( DelaunayTriangulation_src
  src/DelaunayTriangulation.cxx
  src/itkPointSetToDelaunayTriangulationFilter.hxx  
  src/itkPointSetToMeshFilter.hxx
  ${predicates.c_location}
  )

# add header
#
set( DelaunayTriangulation_header
  ${WalkInTriangulation_DIR}/itkWalkInTriangulationFunction.h
  ${PointInCircle_DIR}/itkPointInCircleGeometricalPredicateFunctor.h
  src/itkPointSetToDelaunayTriangulationFilter.h
  src/itkPointSetToMeshFilter.h
  )

# if generate by Xcode add headers to executable
#
if( ${CMAKE_GENERATOR} MATCHES "Xcode" )
  set( DelaunayTriangulation_src
    ${DelaunayTriangulation_src}
    ${DelaunayTriangulation_header}
    )
endif()

# create the executable that will be used by all tests
#
add_executable( DelaunayTriangulation
  ${DelaunayTriangulation_src}
  )

target_link_libraries( DelaunayTriangulation
  ${ITK_LIBRARIES}
  )

# test
#
include( CTest )

add_test( DelaunayTriangulation01
  DelaunayTriangulation regular 5
  )
add_test( DelaunayTriangulation02
  DelaunayTriangulation regular 10
  )
add_test( DelaunayTriangulation03
  DelaunayTriangulation random 10
  )
add_test( DelaunayTriangulation04
  DelaunayTriangulation random 1000
  )