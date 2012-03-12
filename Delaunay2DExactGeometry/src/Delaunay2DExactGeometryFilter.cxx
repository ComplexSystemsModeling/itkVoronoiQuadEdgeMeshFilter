
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"

#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"

#include "itkVTKPolyDataWriter.h"

#include <iostream>
#include <limits>

template< class TMesh >
void CreateDummyTriangle( typename TMesh::Pointer mesh )
{
  typedef TMesh                       MeshType;
  typedef typename MeshType::CellType CellType;

  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;

  if( mesh->GetNumberOfPoints() )
    {
    mesh->Clear();
    mesh->ClearFreePointAndCellIndexesLists();
    }

  int expectedNumPts = 3;
  int expectedNumCells = 1;
  int simpleTriangleCells[3] = 
  { 0, 2, 1 };

  double inf = std::numeric_limits< double >::max();

  typedef typename TMesh::PointType PointType;
  std::vector< PointType > pts( 3 );
  int i( 0 );
  pts[i][0] = -inf; pts[i][1] = -inf; pts[i++][2] = 0.;
  pts[i][0] =  inf; pts[i][1] =  inf; pts[i++][2] = 0.; 
  pts[i][0] =  inf; pts[i][1] = -inf; pts[i++][2] = 0.; 
   
  for( i = 0; i < expectedNumPts; i++ )
    {
    mesh->SetPoint( i, pts[i] );
    }

  typename CellType::CellAutoPointer cellpointer;
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

//
// Delaunay Construction
//

template<
  class TInputPointSet,
  class TOutputMesh
  >
TOutputMesh*
DelaunayTriangulation( TInputPointSet* myPointSet )
{

  typedef          TInputPointSet                          PointSet;
  typedef typename PointSet::PointType                     InputPointType;
  typedef typename PointSet::PointIdentifier               InputPointIdentifier;
  typedef typename PointSet::PointsContainer               InputPointsContainer;
  typedef typename PointSet::PointsContainerConstIterator  InputPointsContainerConstIterator;

  typedef          TOutputMesh                        Mesh;
  typedef typename Mesh::PointType                    PointType;
  typedef typename Mesh::CellType                     CellType;
  typedef typename Mesh::PointIdentifier              PointIdentifier;
  typedef typename Mesh::CellIdentifier               CellIdentifier;
  typedef typename Mesh::PointsContainer              PointsContainer;
  typedef typename Mesh::CellsContainer               CellsContainer;
  typedef typename Mesh::PointIdList                  PointIdList;
  typedef typename Mesh::QEType                       QuadEdgeType;
  typedef typename Mesh::CellsContainerIterator       CellsContainerIteratorType;
  typedef typename Mesh::PointsContainerConstIterator PointsContainerConstIterator;
  
  typedef typename CellType::PointIdConstIterator  PointIdConstIterator;
  typedef typename CellType::PointIdIterator       PointIdIterator;
  typedef typename CellType::CellAutoPointer       CellAutoPointer;

  typedef typename itk::WalkInTriangulationFunction< Mesh >  WalkInTriangulation;
  typedef          itk::VectorContainer< unsigned int, int > CellIdVectorContainer;
  typedef typename itk::QuadEdgeMeshPolygonCell< CellType >  QEPolygonCellType;

  // Create Bounding Triangle
  // Cell ID = 0 and Points ID = {0,1,2}
  typename Mesh::Pointer myMesh = Mesh::New();
  CreateDummyTriangle< Mesh >( myMesh );

  // Loop on all point of the PointSet
  InputPointsContainer *myPoints                  = myPointSet->GetPoints();
  InputPointsContainerConstIterator pointIterator = myPoints->Begin();
  CellIdentifier myStartingCellIndex = 0;

  while( pointIterator != myPoints->End() )
    {
    
    // Get the current point
    InputPointType myTempPoint = pointIterator.Value();
    PointType      myPoint;
    myPoint[0] = myTempPoint[0];
    myPoint[1] = myTempPoint[1];

    // Get the cell id that contain the current point
    CellIdVectorContainer::Pointer cellList = CellIdVectorContainer::New();
    typename WalkInTriangulation::Pointer myWalk = WalkInTriangulation::New();
    try
      {
      cellList = myWalk->Evaluate( myMesh, myPoint, myStartingCellIndex );
      }
    catch( int e )
      {
      std::cerr << "Exception catch : "<< e << std::endl;
      }
    CellIdentifier myCellIndex = cellList->GetElement( cellList->Size()-1 );
    
    // Delete the Cell and replace it by 3 new cell

    // Check Delaunay criterion on the 3 new cell
    
    }
  
  return myMesh;

}


//
// Main function
//

int
main( int argc, char* argv[] )
{
  const unsigned int Dimension = 3;
  typedef double     PixelType;
  
  typedef itk::QuadEdgeMesh< PixelType, Dimension > QEMeshType;
  
  typedef QEMeshType::PointType              PointType;
  typedef QEMeshType::CellType               CellType;
  typedef QEMeshType::CellIdentifier         CellIdentifier; 
  typedef QEMeshType::PointIdentifier        PointIdentifier;
  typedef QEMeshType::CellsContainer         CellsContainer;
  typedef QEMeshType::PointsContainer        PointsContainer;
  typedef QEMeshType::PointIdList            PointIdList;
  typedef QEMeshType::QEType                 QuadEdgeType;
  typedef QEMeshType::CellsContainerIterator CellsContainerIteratorType;

  typedef CellType::PointIdConstIterator     PointIdConstIterator;
  typedef CellType::PointIdIterator          PointIdIterator;
  typedef CellType::CellAutoPointer          CellAutoPointer;

  typedef QuadEdgeType::DualOriginRefType    DualOriginRefType;

  typedef itk::VectorContainer< unsigned int, int > CellVectorContainer;
  typedef itk::WalkInTriangulationFunction< QEMeshType > WalkInTriangulation;

  typedef itk::PointSet< PixelType, Dimension >      PointSetType;
  typedef PointSetType::PointType                    InputPointType;
  typedef PointSetType::PointIdentifier              InputPointIdentifier;
  typedef PointSetType::PointsContainer              InputPointsContainer;
  typedef PointSetType::PointsContainerConstIterator InputPointsContainerIterator;

  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;

  // Input PointSet
  PointSetType::Pointer myPointSet = PointSetType::New();
  InputPointType pts; 
  pts[0] = atof( argv[1] ); 
  pts[1] = atof( argv[2] ); 
  pts[2] = 0;
  myPointSet->SetPoint( 0, pts );

  // New QuadEdgeMesh
  QEMeshType::Pointer myMesh = QEMeshType::New();
  
  CreateDummyTriangle< QEMeshType >( myMesh );

  std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
            << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
            << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
            << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
 
  // Process
  PointType newPts;
  newPts[0] = pts[0]; 
  newPts[1] = pts[1]; 
  newPts[2] = pts[2];
  
  CellIdentifier cellIndex = 0;

  CellVectorContainer::Pointer cellvector = CellVectorContainer::New();
  WalkInTriangulation::Pointer walk = WalkInTriangulation::New();
  try
  {
  cellvector = walk->Evaluate( myMesh, pts, cellIndex );
  }
  catch( int e )
  {
  std::cout << "erreur : " << e << std::endl;
  return EXIT_FAILURE;
  }

  std::cout << "result : "<< cellvector->GetElement( cellvector->Size() -1 ) << std::endl;

  DelaunayTriangulation< PointSetType, QEMeshType >( myPointSet );

  return EXIT_SUCCESS;
}
