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
  #src/Delaunay2DExactGeometryFilter.h
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
  Delaunay2DExactGeometry reg 1
  )
add_test( Delaunay2DExactGeometry02
  Delaunay2DExactGeometry reg 2
  )
add_test( Delaunay2DExactGeometry03
  Delaunay2DExactGeometry reg 3
  )
add_test( Delaunay2DExactGeometry04
  Delaunay2DExactGeometry reg 4
  )
add_test( Delaunay2DExactGeometry05
  Delaunay2DExactGeometry reg 5
  )
add_test( Delaunay2DExactGeometry06
  Delaunay2DExactGeometry reg 6
  )
add_test( Delaunay2DExactGeometry07
  Delaunay2DExactGeometry reg 7
  )
add_test( Delaunay2DExactGeometry08
  Delaunay2DExactGeometry reg 8
  )
add_test( Delaunay2DExactGeometry09
  Delaunay2DExactGeometry reg 9
  )
add_test( Delaunay2DExactGeometry10
  Delaunay2DExactGeometry reg 10
  )

add_test( Delaunay2DExactGeometryRand01
  Delaunay2DExactGeometry rand 20
  )
add_test( Delaunay2DExactGeometryRand02
  Delaunay2DExactGeometry rand 20
  )
add_test( Delaunay2DExactGeometryRand03
  Delaunay2DExactGeometry rand 20
  )
add_test( Delaunay2DExactGeometryRand04
  Delaunay2DExactGeometry rand 20
  )
