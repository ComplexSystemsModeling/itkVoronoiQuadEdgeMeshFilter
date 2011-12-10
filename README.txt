Implementation in progress of a 2D Voronoi/Delaunay algorithm for ITK using Quad Edge Mesh structure.

INPUT: itk::PointSet< TPixelType, VDimension=2 >
OUPUT: itk::QuadEdgeMesh< TPixelType, VDimension=2 >

The algorithm will calculate a Delaunay triangulation and using the dual/primal of the QuadEdgeMesh, the Voronoi tessellation will be generated.

The implementation of the Delaunay will be incremental

# generating a bounding triangle around the input PointSet.
# for all Point p in the PointSet
# find the triangle t that contain p
# 	if p not on an edge
# 		divide t in 3 new triangle t1 t2 t3 such as
# 			t1 = abp, t2 = bcp and t3 = cap
# 		recursively check if Delaunay criterion is respected
# 	else
# 		t' =  neighbor of t that share p
# 		divide t and t' each in 2 new triangle t1, t2, t'1 t'2.
# 		recursively check if Delaunay criterion is respected
# 	endif
# endfor
# remove the bounding triangle and all edge share with it

ToDo :
	Walk in a Triangulation algorithm
	Delaunay criterion recursive check and flip
	Global triangulation algorithm
	Dual/Primal swap

Optional :
	remove/add point from the Delaunay/Voronoi dynamically 
