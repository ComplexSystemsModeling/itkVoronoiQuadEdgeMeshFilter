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
  src/Delaunay2DExactGeometryFilter.cxx
  ${predicates.c_location}
  )

# add header
#
set( Delaunay2DExactGeometry_header
  ${WalkInTriangulation_DIR}/itkWalkInTriangulationFunction.h
  ${PointInCircle_DIR}/itkPointInCircleGeometricalPredicateFunctor.h
  src/Delaunay2DExactGeometryFilter.h
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
  Delaunay2DExactGeometry
  )

