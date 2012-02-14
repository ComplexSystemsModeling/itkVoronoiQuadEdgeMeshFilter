cmake_minimum_required(VERSION 2.8)

set( PointInCircle_DIR "${CMAKE_CURRENT_SOURCE_DIR}/PointInCircle/PointInCircle/" )
set( WalkInTriangulation_DIR "${CMAKE_CURRENT_SOURCE_DIR}/WalkInTriangulation/src" )
#set( DelaunayVoronoi2DExactGeometry_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DelaunayVoronoi2DExactGeometry/" ) # temporary for future development
set( predicates.c_location "${PointInCircle_DIR}/predicates.c" )

enable_testing()
include( CTest )

subdirs( WalkInTriangulation )
#subdirs( DelaunayVoronoi2DExactGeometry ) # temporary for future development
