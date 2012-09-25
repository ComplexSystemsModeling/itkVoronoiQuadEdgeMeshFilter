# How to use the Primal/Dual code
# This small tutorial suposed that you already have CMake and Git installed on your computer.

# 1 - First, get the latest version of ITK from the github repositories of the ComplexSystem.
#     This version contain the latest addition on the itkQuadEdgeMesh structure.
#     Clone git repository

git clone https://github.com/ComplexSystemsModeling/ITK.git

#     Use CMake to configure the project and compile ITK with the Review

# 2 - Second, get the Delaunay Triangulation filter. 
#     Clone the itkVoronoiQuadEdgeMeshFilter repository that contain
#              the itkWalkInTriangulationFunction
#              the itkPointSetToDelaunayTriangulationFilter
#              the 
#     And update the submodule to get the itkPointInCircleFunction

git clone https://github.com/ComplexSystemsModeling/itkVoronoiQuadEdgeMeshFilter.git 
cd itkVoronoiQuadEdgeMeshFilter
git submodule update --init

#     Use CMake to configure and compile the different filters
#     An example of each filter can be found in the source repository of each filter

itkVoronoiQuadEdgeMeshFilter/DelaunayTriangulation/src/DelaunayTriangulation.cxx
itkVoronoiQuadEdgeMeshFilter/WalkInTriangulation/src/WalkInTriangulation.cxx

#     More information on the different filters can be found in the respective Insight Journal

PointInCircle - http://hdl.handle.net/10380/3329
WalkInTriangulation - http://hdl.handle.net/10380/3341
Delaunay Triangulation - http://hdl.handle.net/10380/3372
PrimalToDual - ToBePublished

