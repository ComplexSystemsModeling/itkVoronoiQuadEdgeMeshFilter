cmake_minimum_required( VERSION 2.8 )

project( Delaunay2DExactGeometry )

find_package( ITK REQUIRED) 
include( ${USE_ITK_FILE} ) 

include_directories( ${CMAKE_SOURCE_DIR}/src )
include_directories( ${PointInCircle_DIR} )
include_directories( ${WalkInTriangulation_DIR} )

# add source file
#
set( Delaunay2DExactGeometry_src
  src/DelaunayIncremental.cxx
  src/itkQuadEdgeMeshToDelaunayTriangulationFilter.hxx  
  ${predicates.c_location}
  )

# add header
#
set( Delaunay2DExactGeometry_header
  ${WalkInTriangulation_DIR}/itkWalkInTriangulationFunction.h
  ${PointInCircle_DIR}/itkPointInCircleGeometricalPredicateFunctor.h
  src/itkQuadEdgeMeshToDelaunayTriangulationFilter.h
  )

# if generat by Xcode add headers to executable
#
if( ${CMAKE_GENERATOR} MATCHES "Xcode" )
  set( Delaunay2DExactGeometry_src
    ${Delaunay2DExactGeometry_src}
    ${Delaunay2DExactGeometry_header}
    )
endif()

# create the executable that will be used by all tests
#
add_executable( Delaunay2DExactGeometry
  ${Delaunay2DExactGeometry_src}
  )

target_link_libraries( Delaunay2DExactGeometry
  ${ITK_LIBRARIES}
  )

# test
#
include( CTest )

add_test( Delaunay2DExactGeometry01
  Delaunay2DExactGeometry 1 10
  )
add_test( Delaunay2DExactGeometry02
  Delaunay2DExactGeometry 1 20
  )
add_test( Delaunay2DExactGeometry03
  Delaunay2DExactGeometry 2 50
  )
add_test( Delaunay2DExactGeometry04
  Delaunay2DExactGeometry 2 500
  )
add_test( Delaunay2DExactGeometry05
  Delaunay2DExactGeometry 2 1000
  )
add_test( Delaunay2DExactGeometry06
  Delaunay2DExactGeometry 2 2000
  )
