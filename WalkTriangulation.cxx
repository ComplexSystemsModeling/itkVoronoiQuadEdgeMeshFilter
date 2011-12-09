/*=========================================================================
*
*		Walking in a Triangulation algorithm implementation
*   Based on :
*
*   Devillers, O. and Pion, S. and Teillaud, M. "Walking in a triangulation",
*   Proceedings of the seventeenth annual symposium on Computational geometry,
*   pages 106-114, 2001
*
*   Implementation for ITK by Stéphane Ulysse Rigaud
*		IPAL (Image & Pervasive Access Lab) CNRS - A*STAR
*		Singapore
*
*   Input double  X coordinate of the destination point in the mesh
*         double  Y coordinate of the destination point in the mesh
*
*   Output int I index of the cell that contain the destination point
*          if I = -1, destination point is outside the mesh
*
*   TODO => Include exact discrete geometry predicate
*        => Validation test
*        => Functorise the implementation
*
*=========================================================================*/

#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"

#include "itkVTKPolyDataWriter.h"
#include <iostream>

typedef	double	PixelType;
const	unsigned int	Dimension = 3;

typedef itk::QuadEdgeMesh< PixelType, Dimension >	QEMeshType;

template< class TMesh >
std::vector< typename TMesh::PointType > GeneratePointCoordinates( const unsigned int& iN );

template< class TMesh >
void CreateSquareTriangularMesh( typename TMesh::Pointer mesh );

template< class TPointType >
double orientation ( TPointType r , TPointType a , TPointType b );


int main(int argc, char * argv[] )
{
	if (argc<3) {
	 std::cerr<<"Usage "<<std::endl;
	 std::cerr<<argv[0]<<" DestinationPointX DestinationPointY"<<std::endl;
	 return EXIT_FAILURE;
	}
	
	typedef QEMeshType::PointType	PointType;
	typedef QEMeshType::CellType	CellType;
	
	typedef QEMeshType::PointIdentifier	PointIdentifier;
	typedef QEMeshType::CellIdentifier	CellIdentifier;
	
	typedef QEMeshType::PointsContainer	PointsContainer;
	typedef QEMeshType::CellsContainer	CellsContainer;
	
	typedef QEMeshType::PointIdList	PointIdList;
	typedef QEMeshType::QEType	QuadEdgeType;
	typedef QEMeshType::CellsContainerIterator	CellsContainerIteratorType;
	
	typedef CellType::PointIdConstIterator	PointIdConstIterator;
	typedef CellType::PointIdIterator	PointIdIterator;
	typedef CellType::CellAutoPointer	CellAutoPointer;	
	
	typedef itk::VTKPolyDataWriter< QEMeshType >	MeshWriterType;
	
	QEMeshType::Pointer myMesh = QEMeshType::New();
	CreateSquareTriangularMesh< QEMeshType >	( myMesh );
	
	std::cout<<"\nNo. of Cells : "<<myMesh->GetNumberOfCells()
		<<"\nNo. of Edges : "<<myMesh->GetNumberOfEdges()
		<<"\nNo. of Faces : "<<myMesh->GetNumberOfFaces()
		<<"\nNo. of Points : "<<myMesh->GetNumberOfPoints()<<"\n\n";
	
	MeshWriterType::Pointer writer = MeshWriterType::New();
	
	writer->SetFileName("./myMesh.vtk");
	writer->SetInput(myMesh);
	writer->Update();

	unsigned int orientationTestCompter =0, 
		triangleVisitedCompter =0;
	
	PointType destination;
	destination[0] = atof(argv[1]); destination[1] = atof(argv[2]);
	//destination[0] = 0.75; destination[1] = 3.25;
	
	PointType	pointQ, pointA, pointB, pointC;
	PointIdentifier	pointIdQ, pointIdA, pointIdB, pointIdC;	
	
	CellAutoPointer	myCellPointer;  
	CellsContainer	*myCellsContainer = myMesh->GetCells();
	CellsContainerIteratorType	myCellIterator = myCellsContainer->Begin();
	CellIdentifier	myCellIndex = myCellIterator.Index();
	CellIdentifier	myOldCellIndex;
	
	if (myMesh->GetCell( myCellIndex, myCellPointer)) { 
		triangleVisitedCompter+=1;
		
		PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
		
		pointIdQ = *pointIdIterator;
		myMesh->GetPoint( *pointIdIterator, &pointQ );
		pointIdIterator++;
		
		pointIdB = *pointIdIterator;
		myMesh->GetPoint( *pointIdIterator, &pointB );
		pointIdIterator++;
		
		pointIdC = *pointIdIterator;
		myMesh->GetPoint( *pointIdIterator, &pointC );
		
		if ( orientation( pointB, pointQ, destination ) < 0 )  {
			orientationTestCompter += 1;
			while ( orientation( pointC, pointQ, destination ) < 0 ) {
				orientationTestCompter += 1;
				// r = l
				pointIdB = pointIdC;
				pointB = pointC;
				// t = neighbour( t through ql )
				if (myMesh->FindEdge( pointIdQ, pointIdC )->IsAtBorder()) {
					std::cout<<"This is a border edge, the point is out of the mesh\n";
					myCellIndex = -1; break;
				}
				else {
					triangleVisitedCompter+=1;
					myOldCellIndex = myCellIndex;
					QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetRight();
					QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetLeft();
					if (leftCell == myCellIndex) {
						myCellIndex = rightCell;
					}
					else {
						myCellIndex = leftCell;
					}		
					std::cout<<"We go from cell "<<myOldCellIndex<<" to cell "<<myCellIndex<<" "
						<<"trough the edge "<<pointIdQ<<" - "<<pointIdC<<"\n";
				}
				// l = vertex of t, l!=q, l!=r
				if( myMesh->GetCell( myCellIndex, myCellPointer) ) {
					PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
					pointIdC = *pointIdIterator;
					myMesh->GetPoint( *pointIdIterator, &pointC );
					pointIdIterator++;
					while (pointIdC == pointIdB || pointIdC == pointIdQ) {
						if (pointIdIterator != myCellPointer->PointIdsEnd()) {
							pointIdC = *pointIdIterator;
							myMesh->GetPoint( *pointIdIterator, &pointC );
							pointIdIterator++;
						}
					}
				}
			}
		}
		else {
			orientationTestCompter += 1;
			do {
				orientationTestCompter += 1;
				// l = r
				pointIdC = pointIdB;
				pointC = pointB;
				// t = neighbour( t through ql )
				if (myMesh->FindEdge( pointIdQ, pointIdB )->IsAtBorder()) {
					std::cout<<"This is a border edge, the point is out of the mesh\n";
					myCellIndex = -1; break;
				}
				else {
					triangleVisitedCompter+=1;
					myOldCellIndex = myCellIndex;
					QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetRight();
					QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetLeft();
					if (leftCell == myCellIndex) {
						myCellIndex = rightCell;
					}
					else {
						myCellIndex = leftCell;
					}		
					std::cout<<"We go from cell "<<myOldCellIndex<<" to cell "<<myCellIndex<<" "
						<<"trough the edge "<<pointIdQ<<" - "<<pointIdB<<"\n";
				}
				// r = vertex of t, r!=q, r!=l
				if( myMesh->GetCell( myCellIndex, myCellPointer) ) {
					PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
					pointIdB = *pointIdIterator;
					myMesh->GetPoint( *pointIdIterator, &pointB );
					pointIdIterator++;
					while (pointIdB == pointIdC || pointIdB == pointIdQ) {
						if (pointIdIterator != myCellPointer->PointIdsEnd()) {
							pointIdB = *pointIdIterator;
							myMesh->GetPoint( *pointIdIterator, &pointB );
							pointIdIterator++;
						}
					}
				}
			}
			while ( orientation( pointB, pointQ, destination ) < 0 );			
		}
		// End of initialisation step
		// Q-destination vector has B on its right and C on its left
		while ( orientation( destination, pointB, pointC ) < 0 ) {
			orientationTestCompter += 1;
			// t = neighbour( t through rl )
			if (myMesh->FindEdge( pointIdB, pointIdC )->IsAtBorder()) {
				std::cout<<"This is a border edge, the point is out of the mesh\n";
				myCellIndex = -1; break;
			}
			else {
				triangleVisitedCompter+=1;
				myOldCellIndex = myCellIndex;
				QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdB, pointIdC )->GetRight();
				QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdB, pointIdC )->GetLeft();
				if (leftCell == myCellIndex) {
					myCellIndex = rightCell;
				}
				else {
					myCellIndex = leftCell;
				}
				std::cout<<"We go from cell "<<myOldCellIndex<<" to cell "<<myCellIndex<<" "
					<<"trough the edge "<<pointIdB<<" - "<<pointIdC<<"\n";
			}
			// s = vertex of t, s!=r, s!=l
			if( myMesh->GetCell( myCellIndex, myCellPointer) ) {
				PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
				pointIdA = *pointIdIterator;
				myMesh->GetPoint( *pointIdIterator, &pointA );
				pointIdIterator++;
				while (pointIdA == pointIdB || pointIdA == pointIdC) {
					if (pointIdIterator != myCellPointer->PointIdsEnd()) {
						pointIdA = *pointIdIterator;
						myMesh->GetPoint( *pointIdIterator, &pointA );
						pointIdIterator++;
					}
				}
			}
			if ( orientation( pointA, pointQ, destination ) < 0 ) {
				orientationTestCompter += 1;
				// r = s
				pointIdB = pointIdA;
				pointB = pointA;
			}
			else {
				orientationTestCompter += 1;
				// l = s
				pointIdC = pointIdA;
				pointC = pointA;
			}
		}
		// destination reached
		std::cout<<"We arrived at destination : "<<myCellIndex<<"\n";
		std::cout<< orientationTestCompter <<" orientation test was made \n";
		std::cout<< triangleVisitedCompter <<" triangle was visited \n";
	}
	
	
	return EXIT_SUCCESS;
}

//////////////////////////////////////
// Test Orientation => TODO exact discrete geometry ?

template< class TPointType >
double orientation ( TPointType r , TPointType a , TPointType b )
{	
	double scalar = (a[0] - r[0]) * (b[1] - r[1]) - (a[1] - r[1]) * (b[0] - r[0]);	
	return ( scalar ) >=0 ? 1 : -1 ; 
}

//////////////////////////////////////
// Test Mesh Generation Function

template< class TMesh >
std::vector< typename TMesh::PointType >
GeneratePointCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType        PointType;
  typedef typename PointType::CoordRepType CoordRepType;
  std::vector< PointType > oPt( iN * iN );
	
  for( unsigned int i = 0; i < iN; i++ )
  {
    for( unsigned int j = 0; j < iN; j++ )
    {
      oPt[ i * iN + j ][0] = static_cast< CoordRepType >( j );
      oPt[ i * iN + j ][1] = static_cast< CoordRepType >( i );
      oPt[ i * iN + j ][2] = static_cast< CoordRepType >( 0. );
    }
  }
	
  return oPt;
}

template< class TMesh >
void CreateSquareTriangularMesh( typename TMesh::Pointer mesh )
{
  typedef TMesh                         MeshType;
  typedef typename MeshType::CellType   CellType;
	
  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;
	
  if( mesh->GetNumberOfPoints( ) )
	{
    mesh->Clear( );
    mesh->ClearFreePointAndCellIndexesLists();
	}
	
  /////////////////////////////////////////////////////////////
  int expectedNumPts = 25;
  int expectedNumCells = 32;
  int simpleSquareCells[96] =
  {  0,  1,  6,
		0,  6,  5,
		1,  2,  7,
		1,  7,  6,
		2,  3,  8,
		2,  8,  7,
		3,  4,  9,
		3,  9,  8,
		5,  6, 11,
		5, 11, 10,
		6,  7, 12,
		6, 12, 11,
		7,  8, 13,
		7, 13, 12,
		8,  9, 14,
		8, 14, 13,
    10, 11, 16,
    10, 16, 15,
    11, 12, 17,
    11, 17, 16,
    12, 13, 18,
    12, 18, 17,
    13, 14, 19,
    13, 19, 18,
    15, 16, 21,
    15, 21, 20,
    16, 17, 22,
    16, 22, 21,
    17, 18, 23,
    17, 23, 22,
    18, 19, 24,
    18, 24, 23 };
	
  typedef typename TMesh::PointType PointType;
  std::vector< PointType > pts = GeneratePointCoordinates< TMesh >( 5 );
	
  for(int i=0; i<expectedNumPts; i++)
	{
    mesh->SetPoint( i, pts[i] );
	}
	
  typename CellType::CellAutoPointer cellpointer;
  QEPolygonCellType *poly;
	
  for(int i=0; i<expectedNumCells; i++)
	{
    poly = new QEPolygonCellType( 3 );
    cellpointer.TakeOwnership( poly );
    cellpointer->SetPointId( 0, simpleSquareCells[3*i] );
    cellpointer->SetPointId( 1, simpleSquareCells[3*i+1] );
    cellpointer->SetPointId( 2, simpleSquareCells[3*i+2] );
    mesh->SetCell( i, cellpointer );
	}
}
