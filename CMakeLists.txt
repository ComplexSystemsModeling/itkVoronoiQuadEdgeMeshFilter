cmake_minimum_required(VERSION 2.8.4)

project( Delaunay )

FIND_PACKAGE ( ITK REQUIRED) 
IF ( ITK_FOUND )
  INCLUDE( ${USE_ITK_FILE} ) 
ENDIF( ITK_FOUND )

ADD_EXECUTABLE( Delaunay
  main.cxx
  )
TARGET_LINK_LIBRARIES( Delaunay
${ITK_LIBRARIES}
  )
