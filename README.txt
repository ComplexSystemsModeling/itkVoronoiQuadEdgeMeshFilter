# This work was done by
# Stephane U. Rigaud, Humayun Irshad, Bertrand Moreau and Alexandre Gouaillard
# { rigaud.stephane | humayun.irshad | agouaillard }@gmail.com
# Do not hesitate to contact for any question or problems encountered.
#
# This small tutorial supposed that you are working on a Unix | MacOs | Linux system
# and already have CMake and Git installed on your computer.
# You can find them respectively at

http://www.cmake.org
http://git-scm.com

# 1 - ITK Installation
#     First, open a terminal console and go where you want to have ITK installed.
#     For example

mkdir GitRoot
cd ~/GitRoot/

#     Download the latest version of our modified ITK from the github repositories
#     of the ComplexSystem team. This version contains the latest addition 
#     on the itkQuadEdgeMesh structure that is needed to run our work. The easiest way is to
#     clone the git repository from github where you want to download the sources.

git clone https://github.com/ComplexSystemsModeling/ITK.git

#     A new repository should appeared, called "ITK", that contains all the sources of ITK.
#     You need now to download the possible submodules that are still missing.
#     For that go in the ITK folder and update the submodule

cd ./ITK
git submodule update --init

#     Once the sources and submodule are downloaded, you need to compile ITK.
#     This is totally similar to installing a standard version of ITK.
#     Use CMake and your favorite C++ Compiler.
#     In the ITK folder, create a build directory, and move in it.

mkdir Build
cd ./Build

#     This directory will, in a close future, contain all the compiled library and executable.
#     You need now to generate a makefile that will allow you to compile ITK.
#     For that, CMake is our best friend. Still in the Build repository, run the command
#     ccmake on the CMakeLists.txt located in the ITK directory. In our case it should be

ccmake ..

#     A kind of graphic interface should appear in your terminal that displays a list of 
#     variables and their status.
#     Toggle the "advance mode" by pressing the 't' key (a lot more variable will appear, 
#     keep calm and go on) and set the following variable to ON.

BUILD_EXAMPLES ON
BUILD_TESTING  ON
ITK_USE_REVIEW ON

#     This will assure you that your ITK installation will be totally valid. 
#     However, be prepared to have some card game, coffee, and a good procrastination
#     website as the compilation will take a little bit of time 
#     (intel i5 cpu + 4Gb memory = 2 hours of compilation)  
#     If you have no fear, you can try to win some time by putting OFF the 
#     BUILD_EXAMPLES and the BUILD_TESTING, but we do not advise it.
#
#     It is now time to configure by pressing the 'c' key and then, when the configuration 
#     is done, to press the 'g' key to generate the makefile and leave CMake.
#
#     You should be back to your terminal. Still in the build repository, it is time to 
#     compile everything. For that, type

make

#     And leave him be.
#
#     For more help on the matter (if you are using windows for example), 
#     several tutorial can be found on ITK website

http://www.itk.org/ITK/help/tutorials.html

#     Once ITK is compiled, we should verify that all went well.
#     For that, still in the Build repository, run the command 

ctest

#     This should run automatically a set of test to verify that all is working.
#     It may take a bit of time, but not that much.
#
#     Congratulation! You have successfully installed ITK.
#
# 2 - Delaunay Triangulation and Simplex Mesh project
#     Once ITK is running perfectly, go where you want to install the projects (here the GitRoot folder)
#     and download our sources from github. 

cd ~/GitRoot
git clone https://github.com/ComplexSystemsModeling/itkVoronoiQuadEdgeMeshFilter.git

#     Like previously, the itkVoronoiQuadEdgeMeshFilter repository has been created.
#     Go in, and update possible submodules

cd ./itkVoronoiQuadEdgeMeshFilter
git submodule update --init

#     When is done, you should have in the directory
#       - PointInCircle folder
#       - WalkInTriangulation folder
#       - DelaunayTriangulation folder
#       - SimplexMesh folder
#       - Documentation folder
#       - CMakeLists.txt file
#       - README.txt file
#
#     Like for ITK, create a Build folder

mkdir ./Build
cd ./Build

#     And again, use CMake on the CMakeLists.txt to configure 
#     and compile the different projects.

ccmake ..

#     You should only need to modify the ITK_DIR variable and put the path to the ITK Build 
#     In our case it should be

ITK_DIR ~/GitRoot/ITK/Build

#     Press 'c' and then 'g' key to configure and generate the makefile.
#     And compile

make

#     run the command 

ctest

#     Again, all the tests should pass.
#     Congratulation! You have install our projects on your computer.
#     
#     You should have, in the Build directory, three new directories
#     WalkInTriangulation, DelaunayTriangulation and SimplexMesh
#     That respectfully contain executable of the same name.
#     You can run them to have a quick example of their results.
#
#     An example on how to use each filter can be found in the source repository of 
#     each project.

itkVoronoiQuadEdgeMeshFilter/SimplexMesh/src/SimplexMesh.cxx
itkVoronoiQuadEdgeMeshFilter/DelaunayTriangulation/src/DelaunayTriangulation.cxx
itkVoronoiQuadEdgeMeshFilter/WalkInTriangulation/src/WalkInTriangulation.cxx

#     More information on the different filters can be found 
#     in their respective Insight Journal.

PointInCircle - http://hdl.handle.net/10380/3329
WalkInTriangulation - http://hdl.handle.net/10380/3341
Delaunay Triangulation - http://hdl.handle.net/10380/3372
SimplexMesh - ToBePublished

# Please, be aware that those project are, for some of them, still under progress.
# Therefore you may encounter some bugs, if so, feel free to inform us about them.
# Same if you have any questions about our work.

