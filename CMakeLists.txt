cmake_minimum_required(VERSION 2.8)

set( predicates.c_location "${CMAKE_CURRENT_SOURCE_DIR}/PointInCircle/PointInCircle/predicates.c")

enable_testing()
include( CTest )

subdirs( WalkInTriangulation )
