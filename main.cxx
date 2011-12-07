
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"

#include "itkVTKPolyDataWriter.h"
#include <iostream>

// TODO :
// => Inside Triangle Test
// => Validation Test
// => Functor encapsulation

typedef	double	PixelType;
const	unsigned int	Dimension = 3;
typedef itk::QuadEdgeMesh< PixelType, Dimension >	QEMeshType;

template< class TMesh >
std::vector< typename TMesh::PointType > GeneratePointCoordinates( const unsigned int& iN );

template< class TMesh >
void CreateSquareTriangularMesh( typename TMesh::Pointer mesh );

template< class TPointType >
double orientation ( TPointType p , TPointType q , TPointType r );

template< class TPointType >
bool IsInTriangle ( TPointType p, TPointType a, TPointType b, TPointType c );

int main(int argc, char * argv[] )
{
	
	if (argc<3) {
		std::cerr<<"Usage "<<std::endl;
		std::cerr<<argv[0]<<" Px Py"<<std::endl;
		return EXIT_FAILURE;
	}
	
	//-----------------------------------
	// Typedef

	typedef QEMeshType::PointType	PointType;
	typedef QEMeshType::CellType	CellType;
	typedef QEMeshType::PointIdentifier	PointIdentifier;
	typedef QEMeshType::CellIdentifier	CellIdentifier;
	typedef QEMeshType::PointsContainer	PointsContainer;
	typedef QEMeshType::CellsContainer	CellsContainer;
	typedef QEMeshType::PointIdList	PointIdList;
	typedef QEMeshType::QEType	QEdgeType;

	typedef PointType::VectorType	VectorType;
	typedef CellType::PointIdIterator	PointIdIterator;
	typedef CellType::PointIdConstIterator	PointIdConstIterator;
	
	typedef CellType::CellAutoPointer	CellAutoPointer;	
	typedef QEMeshType::CellsContainerIterator	CellsContainerIteratorType;
  
	typedef itk::VTKPolyDataWriter< QEMeshType >	MeshWriterType;
	
	//-----------------------------------
	// Test Mesh declaration
	
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
	
	//-----------------------------------
	// WALK IN A TRIANGULATION ALGORITHM 
	//
	// INPUT = CellIdentifier
	// OUTPUT = CellIdentifier
	//
	// TEMPORARY INITIALISATION 
	// Get cell t the first cell ID in T  ||   
	// Get cell t the last cell ID processed
	//
	//	Test if p inside t
	//	While seeked point p is not inside cell t
	//	do
	//		Calculate point q the barycentre of cell t
	//		Go through all the edge e-r of the cell t
	//		Until
	//			orientation(pqr) >=0 && orientation(pqe) >=0
	//		If e-r is not border edge
	//			Find cell c that contain the edge e-r and is not cell t
	//			Save c in t
	//			Test if p inside t
	//		else
	//			t = -1
	//			break out while
	//	end while
	//	return t	
	//-----------------------------------
	
	// Test point
	PointType seekPoint;
	seekPoint[0] = atof(argv[1]);
	seekPoint[1] = atof(argv[2]);
	std::cout<<" seek point : "<<seekPoint[0]<<" "<<seekPoint[1]<<"\n";
	
	// compteur
	unsigned int	edgeTest = 0,
								triangleTest = 0,
								edgeCrossed = 0;
	
	// Initialisation
	bool found	= false, 
	     edgeFound	= false;
			
	double direction, directionTest;
		
	PointType A,B,C;
	
	PointType	barycentre;
	PointType	cellPoint, edgePointA, edgePointB, previousEdgePointA, previousEdgePointB;
	PointIdentifier	edgePointIdA, edgePointIdB, previousEdgePointIdA, previousEdgePointIdB;	
	
	CellAutoPointer myCellPointer;  
	CellsContainer  *myCellsContainer = myMesh->GetCells();
	CellsContainerIteratorType myCellIterator = myCellsContainer->Begin();
	CellIdentifier myCellIndex = myCellIterator.Index();
	CellIdentifier myOldCellIndex;
	
	// Initialisation
	// If no starting cell
	// Take the first id cell
	// Verify if pts not inside
	// If not, take the first edge tested correct
	
	// Initial test if the point is in the cell
	// Test if the starting cell exist
	std::cout<<"Starting Id Cell is : "<<myCellIndex<<"\n";
	if (myMesh->GetCell( myCellIndex, myCellPointer)) {
		PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
		
		myMesh->GetPoint( *pointIdIterator, &A );
		pointIdIterator++;
		myMesh->GetPoint( *pointIdIterator, &B );
		pointIdIterator++;		
		myMesh->GetPoint( *pointIdIterator, &C );
		
		myMesh->GetPoint( *pointIdIterator, &C );
		found = IsInTriangle(seekPoint, A, B, C);
		std::cout<<"Is this the correct triangle : "<<found<<"\n";
		
		triangleTest++;
		edgeTest+=3;
		
	}
	
	
	// Do until we found the cell or we reach a border
	while( !found )  {
		
		barycentre[0] = 0;
		barycentre[1] = 0;
		
		// Test if the selected Cell id exist
		if( myMesh->GetCell( myCellIndex, myCellPointer ) ) 
		{
			std::cout<<"We are in Cell Id : "<<myCellIndex<<"\n";
			// Test the type of Cell 0=point 1=line 2=triangle 3=square 4=polygone
			if( myCellPointer->GetType() == 4)													
			{ 
				// Calculate the barycentre of the current cell					
				PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
				while (pointIdIterator != myCellPointer->PointIdsEnd()) {
					myMesh->GetPoint( *pointIdIterator, &cellPoint );
					barycentre[0] = barycentre[0] + cellPoint[0];
					barycentre[1] = barycentre[1] + cellPoint[1];
					pointIdIterator++;
				}
				barycentre[0] = barycentre[0] / myCellPointer->GetNumberOfPoints();
				barycentre[1] = barycentre[1] / myCellPointer->GetNumberOfPoints();
				std::cout<<"The baricentre of the cell "<<myCellIndex<<" is : [ "<<barycentre[0]<<" ; "<<barycentre[1]<<" ]\n";
			}
			else
			{
				std::cout<<"err - not the type of cell expected \n";
			}
		}
		else
		{
			std::cout<<"err - the cell ID was not found \n";
		}

		// Determined the Edge of the Cell that cross the direction we need to go
		// Loop on all point of the cell two by two and test the edge
		PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();		
		int cpt = 0;
		edgeFound = false;
		do {
			myMesh->GetPoint( *pointIdIterator, &edgePointA );
			edgePointIdA = *pointIdIterator;
			pointIdIterator++;
			if (pointIdIterator == myCellPointer->PointIdsEnd()) {
				pointIdIterator = myCellPointer->PointIdsBegin();
			}
			myMesh->GetPoint( *pointIdIterator, &edgePointB );
			edgePointIdB = *pointIdIterator;
			
			std::cout<<"-----------\n";
			std::cout<<"We test the edge : "<<edgePointIdA<<" - "<<edgePointIdB<<"\n";
			
			if (  myMesh->FindEdge( edgePointIdA, edgePointIdB ) != myMesh->FindEdge( previousEdgePointIdA, previousEdgePointIdB ) && 
					  myMesh->FindEdge( edgePointIdA, edgePointIdB ) != myMesh->FindEdge( previousEdgePointIdB, previousEdgePointIdA ) ) 
			{
				
				edgeTest++;
				
				direction = orientation(edgePointA, edgePointB, seekPoint);
				directionTest = orientation(edgePointA, edgePointB, barycentre);
				
				if (direction != directionTest) {
					edgeFound = true;
					std::cout<<"orientation result : "<<direction<<" and "<<directionTest<<"\n";
					std::cout<<"edgefound : "<<edgePointIdA<<" - "<<edgePointIdB<<"\n";
				}
				else {
					std::cout<<"orientation result : "<<direction<<" and "<<directionTest<<"\n";
					std::cout<<"not the good edge\n";
				}
				
			}
			else {
				std::cout<<"we come from this edge, we dont check\n";
			}


			cpt++;
		}
		while ( !edgeFound );
				
		// we have the id and pointer of the selected vertices 
		// => edgePointIdA | *edgePointA and edgePointIdB | *edgePointB
		// we have the id and pointer of the current cell => myCellIndex | myCellPointer
		
		myOldCellIndex = myCellIndex;
		
		// if edge A - B is a border edge then p outside mesh
		// else we look the cells that share the edge
		// and update the id
		if (myMesh->FindEdge( edgePointIdA, edgePointIdB )->IsAtBorder()) {
			std::cout<<"This is a border edge, the point is out of the mesh\n";
			myCellIndex = -1;
			break;
		}
		else { 
			QEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( edgePointIdA, edgePointIdB )->GetLeft();
			QEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( edgePointIdA, edgePointIdB )->GetRight();

			if (leftCell == myCellIndex) {
				myCellIndex = rightCell;
			}
			else {
				myCellIndex = leftCell;
			}
			
			myMesh->GetCell( myCellIndex, myCellPointer);
			previousEdgePointIdA = edgePointIdA;
			previousEdgePointIdB = edgePointIdB;

			std::cout<<"We go from cell "<<myOldCellIndex<<" to cell "<<myCellIndex<<"\n";
			std::cout<<"Trough the edge "<<edgePointIdA<<" - "<<edgePointIdB<<"\n";
			
			edgeCrossed++;
			
			// If Point inside
			PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
			myMesh->GetPoint( *pointIdIterator, &A );
			pointIdIterator++;
			myMesh->GetPoint( *pointIdIterator, &B );
			pointIdIterator++;
			myMesh->GetPoint( *pointIdIterator, &C );
			found = IsInTriangle(seekPoint, A, B, C);
			std::cout<<"Is this the correct triangle : "<<found<<"\n";
			
			triangleTest++;
			edgeTest+=3;
			
		}
		std::cout<<"\n\n";
	}
	std::cout<<"End of the walk \n"; 
	std::cout<<"The point "<<seekPoint[0]<<" "<<seekPoint[1]<<" is int the cell "<<myCellIndex<<"\n";
	std::cout<<"---------------------------------------------\n";
	std::cout<<edgeTest<<" edge test was made\n";
	std::cout<<triangleTest<<" triangle test was made\n";
	std::cout<<edgeCrossed<<" walk was made\n";
	
	// TODO => return cell ID or -1 if outside mesh
	// return myCellIndex;
  return EXIT_SUCCESS;

}

//////////////////////////////////////
// Test Orientation => TODO exact discrete geometry ?

template< class TPointType >
double orientation ( TPointType a , TPointType b , TPointType c )
{
	
	//typedef typename TPointType::VectorType VectorType;
	
	//VectorType pq;
	//pq[0] = q[0] - p[0] ;
	//pq[1] = q[1] - p[1] ;
	
	//VectorType pr;
	//pr[0] = r[0] - p[0] ;
	//pr[1] = r[1] - p[1] ;
	
	//double scalar = pq[0] * pr[0] + pq[1] * pr[1]  ; 
	
	double scalar = (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]);
	
	std::cout<<"scalar product : "<<scalar<<"\n";
	double res = ( scalar ) >=0 ? 1 : -1 ; 
	
	return res;
}

template< class TPointType >
bool IsInTriangle ( TPointType p, TPointType a, TPointType b, TPointType c )
{
	double t1 = ((p[0] - b[0]) * (a[1] - b[1]) - (a[0] - b[0]) * (p[1] - b[1])>=0)?1:0;
	double t2 = ((p[0] - c[0]) * (b[1] - c[1]) - (b[0] - c[0]) * (p[1] - c[1])>=0)?1:0;
	double t3 =	((p[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (p[1] - a[1])>=0)?1:0;
	
	bool results = ( (t1 == t2) && (t2 == t3) );
	
	return results;
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
	
  if( mesh->GetNumberOfPoints( ) ) {
    mesh->Clear( );
    mesh->ClearFreePointAndCellIndexesLists();
	}
	
  /////////////////////////////////////////////////////////////
  int expectedNumPts = 25;
  int expectedNumCells = 32;
  int simpleSquareCells[96] =
  { 0,  1,  6, 
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
	
  for(int i=0; i<expectedNumPts; i++) {
		mesh->SetPoint( i, pts[i] );
	}
	
  typename CellType::CellAutoPointer cellpointer;
  QEPolygonCellType *poly;
	
  for(int i=0; i<expectedNumCells; i++) {
		poly = new QEPolygonCellType( 3 );
    cellpointer.TakeOwnership( poly );
    cellpointer->SetPointId( 0, simpleSquareCells[3*i] );
    cellpointer->SetPointId( 1, simpleSquareCells[3*i+1] );
    cellpointer->SetPointId( 2, simpleSquareCells[3*i+2] );
    mesh->SetCell( i, cellpointer );
	}
}
