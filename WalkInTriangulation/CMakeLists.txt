cmake_minimum_required(VERSION 2.8.4)

project( WalkInTriangulation )

FIND_PACKAGE ( ITK REQUIRED) 
INCLUDE( ${USE_ITK_FILE} ) 

# create the executable that will be used by all tests
ADD_EXECUTABLE( WalkInTriangulation
  src/WalkInTriangulation.cxx
  src/predicates.c
  )
TARGET_LINK_LIBRARIES( WalkInTriangulation
${ITK_LIBRARIES}
  )
