cmake_minimum_required(VERSION 2.8)

set( PointInCircle_DIR "${CMAKE_CURRENT_SOURCE_DIR}/PointInCircle/PointInCircle/" )
set( WalkInTriangulation_DIR "${CMAKE_CURRENT_SOURCE_DIR}/WalkInTriangulation/src" )
set( DelaunayTriangulation_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DelaunayTriangulation/src" )
set( SimplexMesh_DIR "${CMAKE_CURRENT_SOURCE_DIR}/SimplexMesh/src" )
set( predicates.c_location "${PointInCircle_DIR}/predicates.c" )

enable_testing()
include( CTest )

subdirs( WalkInTriangulation )
subdirs( DelaunayTriangulation )
subdirs( SimplexMesh )
