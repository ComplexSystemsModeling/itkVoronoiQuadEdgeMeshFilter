# This work was done by Stéphane U. Rigaud, Humayun Irshad, Bertrand Moreau and Alexandre Gouaillard 
# www.ipal.cnrs.fr

# How to use the Primal/Dual code
# This small tutorial supposed that you already have CMake and Git installed on your computer.
# You can find them respectively at

http://www.cmake.org
http://git-scm.com

# 1 - ITK Installation
#     First, download the latest version of our modified ITK from the github repositories 
#     of the ComplexSystem team. This version contains the latest addition 
#     on the itkQuadEdgeMesh structure that is needed to run our work. The easiest way is to 
#     clone the git repository from github where you want to download the sources.

git clone https://github.com/ComplexSystemsModeling/ITK.git
cd ./ITK
git submodule update --init

#     Once the sources and submodule are downloaded, you need to compile it.
#     This is totally similar to installing a standard version of ITK.
#     Use CMake and your favorite C++ Compiler.
#     For more help on the matter, several tutorial can be found on ITK website

http://www.itk.org/ITK/help/tutorials.html

#     For a better use of ITK, make sure to build examples, testing and reviews 
#     by setting on the following variable in CMake configuration in CMake advance mode.
#     Take note that setting the examples and the testing will make 
#     the compilation much slower, but will unsure that you have a perfectly 
#     installed version of ITK.

BUILD_EXAMPLES ON
BUILD_TESTING  ON
ITK_USE_REVIEW ON

#     Once ITK compiled, run the test command "CTest" into your Build repository
#     in order to check if the installation went correctly.

# 2 - Delaunay Triangulation and Simplex Mesh
#     Once ITK is running perfectly, download our repository from github. 

git clone https://github.com/ComplexSystemsModeling/itkVoronoiQuadEdgeMeshFilter.git 
cd ./itkVoronoiQuadEdgeMeshFilter
git submodule update --init

#     The repository should contain the sources of
#       - the Point In Circle
#       - the Walk In Triangulation
#       - the Delaunay Triangulation
#       - the Simplex Mesh
#       - A CMakeLists.txt
#
#     Use CMake to configure and compile the different projects and 
#     run the command "CTest" into your Build repository to verify that all is good.
#     All the tests should pass.
#     An example of each filter can be found in the source repository of each project.

itkVoronoiQuadEdgeMeshFilter/SimplexMesh/src/SimplexMesh.cxx
itkVoronoiQuadEdgeMeshFilter/DelaunayTriangulation/src/DelaunayTriangulation.cxx
itkVoronoiQuadEdgeMeshFilter/WalkInTriangulation/src/WalkInTriangulation.cxx

#     More information on the different filters can be found 
#     in their respective README file and in their respective Insight Journal.

PointInCircle - http://hdl.handle.net/10380/3329
WalkInTriangulation - http://hdl.handle.net/10380/3341
Delaunay Triangulation - http://hdl.handle.net/10380/3372
PrimalToDual - ToBePublished

# Please, be aware that those project are, for some of them, still under progress.
# Therefore you may encounter some bugs, if so, feel free to inform us about them.

