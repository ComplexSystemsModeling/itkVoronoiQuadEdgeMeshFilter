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
*   IPAL (Image & Pervasive Access Lab) CNRS - A*STAR
*   Singapore
*
*   Input double  X coordinate of the destination point in the mesh
*         double  Y coordinate of the destination point in the mesh
*         unsigned int I starting cell id
*
*   Output int I index of the cell that contain the destination point
*          if I = -1, destination point is outside the mesh
*
*   TODO => Validation test
*        => Functorise the implementation
*
*=========================================================================*/

#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"

#include "vnl/vnl_det.h" 
#include "vnl/vnl_matrix_fixed.h" 

#include "itkVTKPolyDataReader.h"
#include <iostream>

#include "itkWalkInTriangulationFunctor.h"


typedef	double	PixelType;
const	unsigned int	Dimension = 3;

typedef itk::QuadEdgeMesh< PixelType, Dimension >	QEMeshType;

template< class TMesh >
std::vector< typename TMesh::PointType > GeneratePointCoordinates( const unsigned int& iN );

template< class TMesh >
void CreateSquareTriangularMesh( typename TMesh::Pointer mesh );

int main(int argc, char * argv[] )
{
	//if (argc<3) {
	// std::cerr<<"Usage "<<std::endl;
	// std::cerr<<argv[0]<<" DestinationPointX DestinationPointY InitialCellId"<<std::endl;
	// return EXIT_FAILURE;
	//}
	
	// Typedef
	
	typedef QEMeshType::PointType              PointType;
	typedef QEMeshType::CellType               CellType;
	typedef QEMeshType::PointIdentifier        PointIdentifier;
	typedef QEMeshType::CellIdentifier         CellIdentifier;
	typedef QEMeshType::PointsContainer        PointsContainer;
	typedef QEMeshType::CellsContainer         CellsContainer;
	typedef QEMeshType::PointIdList            PointIdList;
	typedef QEMeshType::QEType                 QuadEdgeType;
	typedef QEMeshType::CellsContainerIterator CellsContainerIteratorType;
	
	typedef CellType::PointIdConstIterator     PointIdConstIterator;
	typedef CellType::PointIdIterator          PointIdIterator;
	typedef CellType::CellAutoPointer          CellAutoPointer;	
	
	typedef itk::VTKPolyDataReader< QEMeshType >	MeshReaderType;
	
	// Mesh Initialisation
	
	QEMeshType::Pointer mesh = QEMeshType::New();
	
	MeshReaderType::Pointer reader = MeshReaderType::New();
	reader->SetFileName("StandartTestMesh.vtk");
	reader->Update();
	mesh = reader->GetOutput();
	
	std::cout<<"\nNo. of Cells : "<<mesh->GetNumberOfCells()
	         <<"\nNo. of Edges : "<<mesh->GetNumberOfEdges()
	         <<"\nNo. of Faces : "<<mesh->GetNumberOfFaces()
           <<"\nNo. of Points : "<<mesh->GetNumberOfPoints()<<"\n\n";

	// Walk in a Triangulation Algorithm Test 
	// Standard Mesh Test
	
	PointType pts1, pts2;
	CellIdentifier cell1, cell2;
	CellIdentifier res1, res2;
	
	pts1[0] = 3.75; pts1[1] = 3.25; cell1 = 0;
	pts2[0] = 0.75; pts2[1] = 0.25; cell2 = 31;
	
	res1 = WalkInTriangulation< QEMeshType >( mesh, pts1, cell1, true);
	res2 = WalkInTriangulation< QEMeshType >( mesh, pts2, cell2, true);
	
	// Random Mesh Test
	
	reader->SetFileName("RandomTestMesh.vtk");
	reader->Update();
	mesh = reader->GetOutput();
	
	std::cout<<"\nNo. of Cells : "<<mesh->GetNumberOfCells()
	<<"\nNo. of Edges : "<<mesh->GetNumberOfEdges()
	<<"\nNo. of Faces : "<<mesh->GetNumberOfFaces()
	<<"\nNo. of Points : "<<mesh->GetNumberOfPoints()<<"\n\n";
	
	pts1[0] = 8.75; pts1[1] = 7.25; cell1 = 3;
	pts2[0] = 2.75; pts2[1] = 1.25; cell2 = 10;
	
	res1 = WalkInTriangulation< QEMeshType >( mesh, pts1, cell1, true);
	res2 = WalkInTriangulation< QEMeshType >( mesh, pts2, cell2, true);
	
	return EXIT_SUCCESS;
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
