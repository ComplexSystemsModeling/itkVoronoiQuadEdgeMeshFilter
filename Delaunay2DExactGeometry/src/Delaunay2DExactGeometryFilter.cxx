
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
  typedef          TMesh              MeshType;
  typedef typename MeshType::CellType CellType;
  
  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;

  std::cout<< "dummy triangle creation\n";

  if( mesh->GetNumberOfPoints() )
    {
    mesh->Clear();
    mesh->ClearFreePointAndCellIndexesLists();
    }

  int expectedNumPts = 4;
  int expectedNumCells = 2;
  int simpleTriangleCells[6] = { 0, 1, 3, 0, 3, 2 };

  double max = std::numeric_limits< double >::max(); 
  double min = std::numeric_limits< double >::min();

  min = -10; max = 10;

  typedef typename TMesh::PointType PointType;
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

template<
  class TInputMesh
  >
typename TInputMesh::CellIdentifier
DelaunayRecursiveCriterionEvaluation( 
               TInputMesh* myMesh, 
               typename TInputMesh::PointIdentifier myPointIndex, 
               typename TInputMesh::CellIdentifier myCellIndex )
{

  typedef          TInputMesh                         Mesh;
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

  typedef typename QuadEdgeType::DualOriginRefType DualOriginRefType;

  typedef typename itk::WalkInTriangulationFunction< Mesh >  WalkInTriangulation;
  typedef typename itk::QuadEdgeMeshPolygonCell< CellType >  QEPolygonCellType;

  PointIdentifier pA, pB, pC;
  CellIdentifier t1, t2;
  DualOriginRefType neighbourCellIndex;

  CellAutoPointer myCellPointer;
  myMesh->GetCell( myCellIndex, myCellPointer );
  PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();

  // Get the neighbour cell
  while( pointIdIterator != myCellPointer->PointIdsEnd() )
    {
    pA = *pointIdIterator;
    pB = *(++pointIdIterator);
    if( myPointIndex != pA && myPointIndex != pB )
      {
      if( myMesh->FindEdge( pA, pB )->IsAtBorder() )
        {
        // Cell at border of the mesh
	// No possible check of the delaunay criterion
	std::cout<<"border cell\n";
	return myCellIndex;
	}
      else
        {
        neighbourCellIndex = myMesh->FindEdge( pA, pB )->GetRight();
	}
      }
    }

  // Test of the delaunay criterion using PointInCircle test function
  // from B. Moreau insigh journal 
  PointType myPointCoord = myMesh->GetPoint( myPointIndex );
  if( TestPointInTriangleInMesh< Mesh >( myMesh, neighbourCellIndex, myPointCoord, true ) )
    {
   
    std::cout<<"we flip the edge\n";
	    
    // Criterion fail, need to flip the edge
    // We delete the 2 cells
    myMesh->GetCell( neighbourCellIndex, myCellPointer );
    pointIdIterator = myCellPointer->PointIdsBegin();
    while( pointIdIterator != myCellPointer->PointIdsEnd() )
      {
      if( pA != *pointIdIterator && pB != *pointIdIterator )
        {
        pC = *pointIdIterator;
	}
      ++pointIdIterator;
      }
    myMesh->DeleteFace( myCellIndex );
    myMesh->DeleteFace( neighbourCellIndex );

    // We create 2 new cells
    int simpleTriangleCells[6] = 
    { myPointIndex, pC, pA, 
      myPointIndex, pB, pC };
    
    QEPolygonCellType *poly;
    CellAutoPointer cellpointer;
    CellIdentifier cellIndexTab[2];
    for( unsigned int i = 0; i<2; i++ )
      {
      poly = new QEPolygonCellType( 3 );
      cellpointer.TakeOwnership( poly );
      cellpointer->SetPointId( 0, simpleTriangleCells[i*3] );
      cellpointer->SetPointId( 1, simpleTriangleCells[i*3+1] );
      cellpointer->SetPointId( 2, simpleTriangleCells[i*3+2] );
      cellIndexTab[i] = myMesh->FindFirstUnusedCellIndex();
      myMesh->SetCell( cellIndexTab[i], cellpointer ); 
      } 
    
    // we check if the 2 new cells respect the criterion
    myCellIndex = DelaunayRecursiveCriterionEvaluation( myMesh, myPointIndex, cellIndexTab[0] ); 
    DelaunayRecursiveCriterionEvaluation( myMesh, myPointIndex, cellIndexTab[1] ); 
    }

  return myCellIndex;
}

//
// Delaunay Construction
//

template<
  class TInputPointSet,
  class TOutputMesh
  >
typename TOutputMesh::Pointer
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

  typedef itk::QuadEdgeMesh< double, 3> TTMesh;

  // Create Bounding Triangle
  // Cell ID = 0 and Points ID = {0,1,2}
  typename Mesh::Pointer myMesh = Mesh::New();
  CreateDummyTriangle< Mesh >( myMesh );
  std::cout << "\nStarting Mesh\n";
  std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
            << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
            << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
            << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n"; 

  // Loop on all point of the PointSet
  InputPointsContainer *myPoints                  = myPointSet->GetPoints();
  InputPointsContainerConstIterator pointIterator = myPoints->Begin();
  CellIdentifier myStartingCellIndex = 0; // get the first cell instead of 0

  int cpt = 0;
  while( pointIterator != myPoints->End() )
    {
     
    // Get the current point
    InputPointType myTempPoint = pointIterator.Value();
    PointType      myPoint;
    myPoint[0] = myTempPoint[0];
    myPoint[1] = myTempPoint[1];
    myPoint[2] = 0.;
    PointIdentifier myPointIndex = myMesh->FindFirstUnusedPointIndex();

    std::cout << "--------------------------------------"<< std::endl;
    std::cout << "Iteration point : "<< cpt++ 
	      << " - coord: " << myPoint[0] << ";" << myPoint[1] << std::endl;

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
      if( e == -1 )
        { 
        std::cout<< "Error - Point is out of the mesh\n";
        }
      else if( e == -2 )
        {
        std::cerr<< "Error - Starting cell does not exist\n";
        }
      else
        {
	std::cerr<< "Error - Unknown error\n";
	}
      }

    CellIdentifier myCellIndex;
    if( cellList->Size() > 0 )
      {
      myCellIndex = ( --cellList->End() )->Value();
      std::cout<< "Retrived Cell : "<< myCellIndex << std::endl;
      }
    else
      {
      myCellIndex = std::numeric_limits<unsigned int>::max();
      std::cout<< "Error - Retrived Cell : "<< myCellIndex << std::endl;
      }

    // Delete the Cell and replace it by 3 new cell
    CellAutoPointer myCellPointer;
    if( myMesh->GetCell( myCellIndex, myCellPointer ) )
      {
      myMesh->GetCell( myCellIndex, myCellPointer );
      PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
      PointIdentifier pA = *pointIdIterator;
      pointIdIterator++;
      PointIdentifier pB = *pointIdIterator;
      pointIdIterator++;
      PointIdentifier pC = *pointIdIterator;
   
      myMesh->DeleteFace( myCellIndex );
      myMesh->SetPoint( myPointIndex, myPoint );

      CellIdentifier t1 = myMesh->FindFirstUnusedCellIndex();
      myMesh->AddFaceTriangle( pA, pB, myPointIndex );
      std::cout<<"cell "<<t1<<" = "<<pA<<" "<<pB<<" "<<myPointIndex<<std::endl;
      CellIdentifier t2 = myMesh->FindFirstUnusedCellIndex();
      myMesh->AddFaceTriangle( pB, pC, myPointIndex );
      std::cout<<"cell "<<t2<<" = "<<pB<<" "<<pC<<" "<<myPointIndex<<std::endl;
      CellIdentifier t3 = myMesh->FindFirstUnusedCellIndex();
      myMesh->AddFaceTriangle( pC, pA, myPointIndex );
      std::cout<<"cell "<<t3<<" = "<<pC<<" "<<pA<<" "<<myPointIndex<<std::endl;

      // Check Delaunay criterion on the 3 new cell
      //t1 = DelaunayRecursiveCriterionEvaluation< Mesh >( myMesh, myPointIndex, t1 );
      //t2 = DelaunayRecursiveCriterionEvaluation< Mesh >( myMesh, myPointIndex, t2 );
      //t3 = DelaunayRecursiveCriterionEvaluation< Mesh >( myMesh, myPointIndex, t3 );

      std::cout<< "new cell id : "<< t1<<" ; "<<t2<<" ; "<<t3<<"\n";
      std::cout << "\nMesh Update\n";
      std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
                << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
                << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
                << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
    
      // Nest iteration setup
      myStartingCellIndex = t1;
      ++pointIterator;

      //debug
      // 

      typedef typename itk::VTKPolyDataWriter<TTMesh> MeshWriter;
      MeshWriter::Pointer write = MeshWriter::New();
      write->SetFileName("./tempMesh.vtk");
      write->SetInput( myMesh );
      write->Update();  

      getchar();
      }
    else
      {
      std::cout<<"Error - Could not find the specified cell\n";
      break;
      }
    }

  std::cout<< "end of delaunay  " << std::endl;  
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
  
  typedef itk::PointSet< PixelType, Dimension >     PointSetType;
  typedef PointSetType::PointType                   PointType;

  int n = atoi( argv[1] );

  // Input PointSet
  PointSetType::Pointer myPointSet = PointSetType::New();
  std::vector<PointType> pts = GeneratePointCoordinates< PointSetType >( n ); 
  for( int i = 0; i < (n*n) ; i++ )
    {
    myPointSet->SetPoint( i, pts[i] );
    }

  // Dummy QuadEdgeMesh test
  QEMeshType::Pointer myMesh = QEMeshType::New();

  // Delaunay Loop Test
  myMesh = DelaunayTriangulation< PointSetType, QEMeshType >( myPointSet );

  

  return EXIT_SUCCESS;
}
