cmake_minimum_required(VERSION 2.8.4)

project( walkTriangulation )

FIND_PACKAGE ( ITK REQUIRED) 
IF ( ITK_FOUND )
  INCLUDE( ${USE_ITK_FILE} ) 
ENDIF( ITK_FOUND )

ADD_EXECUTABLE( walkTriangulation
  WalkTriangulation.cxx
  )
TARGET_LINK_LIBRARIES( walkTriangulation
${ITK_LIBRARIES}
  )
