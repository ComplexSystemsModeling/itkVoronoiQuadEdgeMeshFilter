
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"

#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"

#include <iostream>
#include <limits>

//TEMPORARY
#include "itkVTKPolyDataWriter.h"

// Create a dummy mesh for triangulation initialisation
// Mesh points should be at inifinity
//
template< class TMeshType >
void CreateDummyMesh( typename TMeshType::Pointer mesh, typename TMeshType::PixelType limit )
{
 // Generate a square triangulated mesh
 //
 //   2---------3
 //   |       / |
 //   |  1  /   |
 //   |   /     |
 //   | /    0  |
 //   0---------1
 //
 // Anti-ClockWise orientation
 // Cell 0 -> 013
 // Cell 1 -> 032
 
 typedef          TMeshType           MeshType;
 typedef typename MeshType::CellType  CellType;
 typedef typename MeshType::PointType PointType;
 typedef typename MeshType::PixelType PixelType;
 
 typedef typename CellType::CellAutoPointer CellAutoPointer;
 
 typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;
 
 if( mesh->GetNumberOfPoints() )
 {
	mesh->Clear();
	mesh->ClearFreePointAndCellIndexesLists();
 }
 
 int expectedNumPts = 4;
 int expectedNumCells = 2;
 int simpleTriangleCells[6] = { 0, 1, 3, 0, 3, 2 };
 
 PixelType min = -limit; // Change to infinity if possible 
 PixelType max =  limit; // Change to infinity if possible
 
 std::vector< PointType > pts( 4 );
 int i( 0 );
 
 pts[i][0] = min; pts[i][1] = min; pts[i++][2] = 0.;
 pts[i][0] = max; pts[i][1] = min; pts[i++][2] = 0.; 
 pts[i][0] = min; pts[i][1] = max; pts[i++][2] = 0.; 
 pts[i][0] = max; pts[i][1] = max; pts[i++][2] = 0.;
 
 for( i = 0; i < expectedNumPts; i++ )
 {
	mesh->SetPoint( i, pts[i] );
 }
 
 CellAutoPointer cellpointer;
 QEPolygonCellType *poly;
 
 for( i = 0; i < expectedNumCells; i++ )
 {
	poly = new QEPolygonCellType( 3 );
	cellpointer.TakeOwnership( poly );
	cellpointer->SetPointId( 0, simpleTriangleCells[i*3] );
	cellpointer->SetPointId( 1, simpleTriangleCells[i*3+1] );
	cellpointer->SetPointId( 2, simpleTriangleCells[i*3+2] );
	mesh->SetCell( i, cellpointer );
 } 
}


// Recursively check delaunay criterion
//
template< class TMeshType >
typename TMeshType::CellIdentifier
DelaunayRecursiveCriterionEvaluation( 
																		 typename TMeshType::Pointer mesh,
																		 typename TMeshType::PointIdentifier myPointIndex, 
																		 typename TMeshType::CellIdentifier myCellIndex 
																		 )
{
 typedef          TMeshType												       MeshType;
 typedef typename MeshType::Pointer                      MeshTypePointer;
 typedef typename MeshType::PointType                    MeshPointType;
 typedef typename MeshType::CellType                     MeshCellType;
 typedef typename MeshType::PointIdentifier              MeshPointIdentifier;
 typedef typename MeshType::CellIdentifier               MeshCellIdentifier;
 typedef typename MeshType::PointsContainer              MeshPointsContainer;
 typedef typename MeshType::CellsContainer               MeshCellsContainer;
 typedef typename MeshType::PointIdList                  MeshPointIdList;
 typedef typename MeshType::QEType                       MeshQuadEdgeType;
 typedef typename MeshType::CellsContainerIterator       MeshCellsContainerIteratorType;
 typedef typename MeshType::PointsContainerConstIterator MeshPointsContainerConstIterator;
 
 typedef typename MeshCellType::PointIdConstIterator  MeshCellPointIdConstIterator;
 typedef typename MeshCellType::PointIdIterator       MeshCellPointIdIterator;
 typedef typename MeshCellType::CellAutoPointer       MeshCellCellAutoPointer;

 typedef typename MeshQuadEdgeType::DualOriginRefType MeshDualOriginRefType;
 
 typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType > QEPolygonCellType;

 MeshPointIdentifier     pA,pB,pC;
 MeshPointIdentifier     p[3];
 MeshPointType           myPointCoord;
 MeshCellIdentifier      t1, t2;
 MeshDualOriginRefType   adjacentCellIndex;
 MeshCellCellAutoPointer myCellPointer;
 MeshCellPointIdIterator pointIdIterator; 
 
 mesh->GetCell( myCellIndex, myCellPointer );
 pointIdIterator = myCellPointer->PointIdsBegin();
 p[0] = *(pointIdIterator);
 p[1] = *(++pointIdIterator);
 p[2] = *(++pointIdIterator);
 
 for( int i = 0; i < 3; i++ )
 {
	if( p[i%3] == myPointIndex )
	{
	 pA = p[(i+1)%3];
	 pB = p[(i+2)%3];	  
	 if( mesh->FindEdge( pA, pB )->IsAtBorder() )
	 {
		return myCellIndex;
	 }
	 else
	 {
		adjacentCellIndex = mesh->FindEdge( pA, pB )->GetRight();
	 }
	}
 }
 
 myPointCoord = mesh->GetPoint( myPointIndex );
 if( TestPointInTriangleInMesh< MeshType >( mesh, adjacentCellIndex, myPointCoord, true ) )
 {

	mesh->GetCell( adjacentCellIndex, myCellPointer );
	pointIdIterator = myCellPointer->PointIdsBegin();
	
	while( pointIdIterator != myCellPointer->PointIdsEnd() )
	{
	 if( pA != *pointIdIterator && pB != *pointIdIterator )
	 {
		pC = *pointIdIterator;
	 }
	 ++pointIdIterator;
	}
	
	mesh->DeleteFace( myCellIndex );
	mesh->DeleteFace( adjacentCellIndex );
	
	int simpleTriangleCells[6] =  // order change the error (segfault/outofmesh)
	{ pA, pC, myPointIndex, 
    pC, pB, myPointIndex };
	
	QEPolygonCellType *poly;
	MeshCellCellAutoPointer cellpointer;
	MeshCellIdentifier      cellIndexTab[2];
	for( unsigned int i = 0; i<2; i++ )
	{
	 poly = new QEPolygonCellType( 3 );
	 cellpointer.TakeOwnership( poly );
	 cellpointer->SetPointId( 0, simpleTriangleCells[i*3] );
	 cellpointer->SetPointId( 1, simpleTriangleCells[i*3+1] );
	 cellpointer->SetPointId( 2, simpleTriangleCells[i*3+2] );
	 cellIndexTab[i] = mesh->FindFirstUnusedCellIndex();
	 mesh->SetCell( cellIndexTab[i], cellpointer ); 
	} 
	
	std::cout<<"\tedge ("<<pA<<"-"<<pB<<") is flip, we delete cells "<<myCellIndex<<" and "<<adjacentCellIndex<<"\n"
	         <<"\t\t new cell : "<<cellIndexTab[0]<<" ("<<myPointIndex<<","<<pA<<","<<pC<<")\n"
	         <<"\t\t new cell : "<<cellIndexTab[1]<<" ("<<myPointIndex<<","<<pC<<","<<pB<<")\n";
		
	//DelaunayRecursiveCriterionEvaluation< MeshType >( mesh, myPointIndex, cellIndexTab[0] ); 
	//DelaunayRecursiveCriterionEvaluation< MeshType >( mesh, myPointIndex, cellIndexTab[1] );
	
	myCellIndex = cellIndexTab[0];
 }

 return myCellIndex;
}

// Delaunay triangulation process loop
//
template<
class TPointSetType,
class TMeshType
>
typename TMeshType::Pointer
DelaunayTriangulation( TPointSetType* myPointSet )
{
 typedef          TPointSetType                               PointSetType;
 typedef typename PointSetType::PointType                     PointSetPointType;
 typedef typename PointSetType::PointIdentifier               PointSetPointIdentifier;
 typedef typename PointSetType::PointsContainer               PointSetPointsContainer;
 typedef typename PointSetType::PointsContainerConstIterator  PointSetPointsContainerConstIterator;
 
 typedef          TMeshType												       MeshType;
 typedef typename MeshType::Pointer                      MeshTypePointer;
 typedef typename MeshType::PointType                    MeshPointType;
 typedef typename MeshType::CellType                     MeshCellType;
 typedef typename MeshType::PointIdentifier              MeshPointIdentifier;
 typedef typename MeshType::CellIdentifier               MeshCellIdentifier;
 typedef typename MeshType::PointsContainer              MeshPointsContainer;
 typedef typename MeshType::CellsContainer               MeshCellsContainer;
 typedef typename MeshType::PointIdList                  MeshPointIdList;
 typedef typename MeshType::QEType                       MeshQuadEdgeType;
 typedef typename MeshType::CellsContainerIterator       MeshCellsContainerIteratorType;
 typedef typename MeshType::PointsContainerConstIterator MeshPointsContainerConstIterator;
 
 typedef typename MeshCellType::PointIdConstIterator  MeshCellPointIdConstIterator;
 typedef typename MeshCellType::PointIdIterator       MeshCellPointIdIterator;
 typedef typename MeshCellType::CellAutoPointer       MeshCellCellAutoPointer;
 
 typedef typename itk::WalkInTriangulationFunction< MeshType > WalkInTriangulation;
 typedef typename WalkInTriangulation::Pointer                 WalkInTriangulationPointer;
 typedef          itk::VectorContainer< unsigned int, int >    MeshCellIdVectorContainer;
 typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType > QEPolygonCellType;
 
 // Initialisation
 
 MeshTypePointer myMesh = MeshType::New();
 CreateDummyMesh< MeshType >( myMesh, 5000 );
 
 std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
           << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
					 << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
           << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
 
 PointSetPointsContainer              *myPoints           = myPointSet->GetPoints();
 PointSetPointsContainerConstIterator pointIterator       = myPoints->Begin();
 MeshCellIdentifier                   myStartingCellIndex = 0;
 
 // Process Loop
 int o(1);
 while( pointIterator != myPoints->End() ) 
 {
	std::cout<<"/---------------------------------------------------/\n";
	std::cout<<"Delaunay iteration : "<<o<<"/"<<myPoints->Size();
	std::cout<<" - pts "<<o+3<<" ("<<pointIterator.Value()[0]<<","<<pointIterator.Value()[1]<<")"<<std::endl;
	
	PointSetPointType       myTempPoint;
	MeshPointType           myPoint;
	MeshPointIdentifier     myPointIndex;
	MeshCellIdentifier      myCellIndex;
	MeshCellCellAutoPointer myCellPointer;

	// Get the current point
	myTempPoint  = pointIterator.Value();
	myPoint[0]   = myTempPoint[0];
	myPoint[1]   = myTempPoint[1];
	myPoint[2]   = myTempPoint[2];
	myPointIndex = myMesh->FindFirstUnusedPointIndex();
	
	// Find the triangle that include the point
	std::cout<<"I walk the mesh -  starting cell "<<myStartingCellIndex<<"\n";
	WalkInTriangulationPointer         letsWalk       = WalkInTriangulation::New();
	MeshCellIdVectorContainer::Pointer walkCellIdList = MeshCellIdVectorContainer::New();
	try 
	{
	 walkCellIdList = letsWalk->Evaluate( myMesh, myPoint, myStartingCellIndex );
	}
	catch( int e ) 
	{
	 if( e == -1 )
	 { 
		std::cout<< "Error Walk - Point is out of the mesh\n";
	 }
	 else if( e == -2 )
	 {
		std::cerr<< "Error Walk - Starting cell does not exist\n";
	 }
	 else
	 {
		std::cerr<< "Error Walk - Unknown error\n";
	 }	 
	 break;
	}
	myCellIndex = ( --walkCellIdList->End() )->Value();
	std::cout<<"I walked the mesh\n";
	
	if( myMesh->GetCell( myCellIndex, myCellPointer ) )
	{	 
	 MeshCellPointIdIterator pointIdIterator;
	 MeshPointIdentifier pA, pB, pC;
	 MeshCellIdentifier  t1, t2, t3;
	 
	 myMesh->GetCell( myCellIndex, myCellPointer );
	 
   pointIdIterator = myCellPointer->PointIdsBegin();
   pA = *pointIdIterator;
	 pointIdIterator++;
   pB = *pointIdIterator;
	 pointIdIterator++;
   pC = *pointIdIterator;
   
	 // Delete the triangle
	 myMesh->DeleteFace( myCellIndex );
	 myMesh->SetPoint( myPointIndex, myPoint );

	 // Create 3 new triangles
	 t1 = myMesh->FindFirstUnusedCellIndex();
	 myMesh->AddFaceTriangle( pA, pB, myPointIndex );
   t2 = myMesh->FindFirstUnusedCellIndex();
	 myMesh->AddFaceTriangle( pB, pC, myPointIndex );
   t3 = myMesh->FindFirstUnusedCellIndex();
	 myMesh->AddFaceTriangle( pC, pA, myPointIndex );
	 
	 std::cout<<"Cell id "<<myCellIndex<<" ( "<<pA<<","<<pB<<","<<pC<<") replaced by \n"
	          <<"\t"<<t1<<" ("<<pA<<","<<pB<<","<<myPointIndex<<")\n"
	          <<"\t"<<t2<<" ("<<pB<<","<<pC<<","<<myPointIndex<<")\n"
	          <<"\t"<<t3<<" ("<<pC<<","<<pA<<","<<myPointIndex<<")\n";
	 
	 // Delaunay criterion recursive test
	 t1 = DelaunayRecursiveCriterionEvaluation< MeshType >( myMesh, myPointIndex, t1 );
	 t2 = DelaunayRecursiveCriterionEvaluation< MeshType >( myMesh, myPointIndex, t2 );
	 t3 = DelaunayRecursiveCriterionEvaluation< MeshType >( myMesh, myPointIndex, t3 );
	 
	 
	 //myStartingCellIndex = 0; // => error pts out of mesh inside WiT
	 myStartingCellIndex = t1; // => error segmentation fault at FindEdge inside WiT
	 //myStartingCellIndex = t2; // => Against all odds it is working
   //myStartingCellIndex = t3; // => error pts out of mesh inside WiT
	 
	 ++pointIterator;
	 
	 // VERIFICATION /////////////////////////////////////////////
	 typedef itk::QuadEdgeMesh< double, 3 > TempMeshType;
	 typedef typename itk::VTKPolyDataWriter< TempMeshType > MeshWriter;
	 MeshWriter::Pointer write = MeshWriter::New();
	 write->SetFileName("./tempMesh.vtk");
	 write->SetInput( myMesh );
	 write->Update();
	 //getchar();
	 // VERIFICATION /////////////////////////////////////////////
	}
	else
	{
	 std::cout<<"Error DelaunayLoop - Could not find the specified cell\n";
	 break;
	}  
	o++;
	std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
						<< "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
	          << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
	          << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
	
	// Detele dummy points
	
 } 
 return myMesh;
}

// Main function
//
int
main( int argc, char* argv[] )
{
 const unsigned int Dimension = 3;
 typedef double     PixelType;
 
 typedef itk::QuadEdgeMesh< PixelType, Dimension > MeshType; 
 typedef itk::PointSet< PixelType, Dimension >     PointSetType;
 typedef PointSetType::PointType                   PointType;
 
 // -------------------------------------------------- //
 // Dummy QuadEdgeMesh test
 MeshType::Pointer myDummyMesh = MeshType::New();
 CreateDummyMesh< MeshType >( myDummyMesh, 5000 );

 // -------------------------------------------------- //
 // Delaunay Loop Test
 int meshSize = atoi( argv[1] );
 int expectedNumPts = meshSize*meshSize;
 
 PointSetType::Pointer myPointSet = PointSetType::New();
 std::vector<PointType> pts = GeneratePointCoordinates< PointSetType >( meshSize ); 
 for( int i = 0; i < expectedNumPts ; i++ )
 {
	myPointSet->SetPoint( i, pts[i] );
 }
 
 MeshType::Pointer myTriangulatedMesh = MeshType::New();
 myTriangulatedMesh = DelaunayTriangulation< PointSetType, MeshType >( myPointSet );
 
 return EXIT_SUCCESS;
}
