

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// O ring check function
// TOBEDELETE
//
template< class TMeshType >
void CheckONextLink( typename TMeshType::Pointer mesh )
{
  typedef          TMeshType                         MeshType;
  typedef typename MeshType::PointType               PointType;
  typedef typename MeshType::PointIdentifier         PointIndex;
  typedef typename MeshType::PointsContainer         PointsContainer;
  typedef typename MeshType::PointsContainerIterator PointsContainerIterator;
  typedef typename MeshType::QEType                  QuadEdgeType;
  
  typedef typename itk::VTKPolyDataWriter< MeshType > MeshWriterType;
  typedef typename MeshWriterType::Pointer                MeshWriterPointer;
  
  PointsContainer*        myPoints;
  PointsContainerIterator myPointIterator;
  
  myPoints = mesh->GetPoints();
  myPointIterator = myPoints->Begin();
  
  while( myPointIterator != myPoints->End() )
    {
    PointType p = myPointIterator.Value();
    QuadEdgeType* e = p.GetEdge();
    QuadEdgeType* next = e->GetOnext();
    while( e != next )
      {
      next = next->GetOnext();
      } 
    ++myPointIterator;
    }
  std::cout<<"Check O ring done\n";
}

// Random coordonates generation function
//
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateRandomPointCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType        PointType;
  typedef typename PointType::CoordRepType CoordRepType;
  std::vector< PointType > oPt( iN * iN );
  
  srand(time(NULL));
  for( unsigned int i = 0; i < iN; i++ )
    {
    oPt[ i ][0] = static_cast< CoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][1] = static_cast< CoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][2] = static_cast< CoordRepType >( 0. );
    }
  return oPt;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Dummy Mesh creation function
// Create a dummy square mesh for triangulation initialisation
// limit argument is to determine the size of the square
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
  
  PixelType min = -limit;  
  PixelType max =  limit;
  
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
  
  CheckONextLink<MeshType>( mesh );
}

// Delaunay triangulation process loop
//
template<
class TPointSetType,
class TMeshType
>
typename TMeshType::Pointer
DelaunayTriangulation( typename TPointSetType::Pointer myPointSet )
{
  typedef          TPointSetType                               PointSetType;
  typedef typename PointSetType::PointType                     PointSetPointType;
  typedef typename PointSetType::PointIdentifier               PointSetPointIdentifier;
  typedef typename PointSetType::PointsContainer               PointSetPointsContainer;
  typedef typename PointSetType::PointsContainerConstIterator  PointSetPointsContainerConstIterator;
  
  typedef          TMeshType                              MeshType;
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
  
  CheckONextLink<MeshType>( myMesh );
  
  PointSetPointsContainer              *myPoints           = myPointSet->GetPoints();
  PointSetPointsContainerConstIterator pointIterator       = myPoints->Begin();
  MeshCellIdentifier                   myStartingCellIndex = 0;
  
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
        std::cerr<< "Error Walk - Point is out of the mesh\n";
        throw -1;
        }
      else if( e == -2 )
        {
        std::cerr<< "Error Walk - Starting cell does not exist\n";
        throw -1;
        }
      else
        {
        std::cerr<< "Error Walk - Unknown error\n";
        throw -2;
        }    
      break;
      }
    myCellIndex = ( --walkCellIdList->End() )->Value();
    std::cout<<"I walked the mesh: "<<myCellIndex<<"\n";
    
    if( myMesh->GetCell( myCellIndex, myCellPointer ) )
      { 
      MeshCellPointIdIterator pointIdIterator;
      std::vector< MeshCellIdentifier > cellPointsIds( 3 );
      std::vector< MeshCellIdentifier > newCellIds( 3 );
      MeshCellCellAutoPointer cellpointer;
      QEPolygonCellType *poly;
      
      myMesh->GetCell( myCellIndex, myCellPointer );
      
      pointIdIterator = myCellPointer->PointIdsBegin();
      cellPointsIds[0] = *pointIdIterator;
      pointIdIterator++;
      cellPointsIds[1] = *pointIdIterator;
      pointIdIterator++;
      cellPointsIds[2] = *pointIdIterator;
        
      myMesh->DeleteFace( myCellIndex );     

      myPointIndex = myMesh->FindFirstUnusedPointIndex();
      myMesh->SetPoint( myPointIndex, myPoint );

      // Create 3 new anticlockwise oriented triangles
      for( unsigned int i = 0; i < 3; i++ )
        {
        newCellIds[i] = myMesh->FindFirstUnusedCellIndex();
        poly = new QEPolygonCellType( 3 );
        cellpointer.TakeOwnership( poly );
        cellpointer->SetPointId( 0, cellPointsIds[ (i)   % 3 ] );
        cellpointer->SetPointId( 1, cellPointsIds[ (i+1) % 3 ] );
        cellpointer->SetPointId( 2, myPointIndex );
        myMesh->SetCell( newCellIds[i], cellpointer );
        } 
        
      std::cout << "Cell id " << myCellIndex <<" ( "<< cellPointsIds[0] << "," << cellPointsIds[1] << "," << cellPointsIds[2]; 
      std::cout << ") replaced by \n";
      std::cout << "\t" << newCellIds[0] << " (" << cellPointsIds[0] << "," << cellPointsIds[1] << "," << myPointIndex << ")\n";
      std::cout << "\t" << newCellIds[1] << " (" << cellPointsIds[1] << "," << cellPointsIds[2] << "," << myPointIndex << ")\n";
      std::cout << "\t" << newCellIds[2] << " (" << cellPointsIds[2] << "," << cellPointsIds[0] << "," << myPointIndex << ")\n";  
      }
      
    std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
              << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
              << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
              << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
      
    o++;
    ++pointIterator;
    }    
  
  CheckONextLink<MeshType>( myMesh );
  
  return myMesh;
}


// Main function
//
int
main( int argc, char* argv[] )
{
  if( argc != 3 )
    {
    std::cerr << "Exit wrong arguments\n";
    std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]\n";
    return EXIT_FAILURE;
    }
  
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

  int type = atoi( argv[1] );
  int meshSize = atoi( argv[2] );
  int expectedNumPts = 0;
  std::vector<PointType> pts;

  switch(type) 
  {
    case 1 :     
      pts = GeneratePointCoordinates< PointSetType >( meshSize ); 
      expectedNumPts = meshSize*meshSize;
      break;
    case 2 :
      pts = GenerateRandomPointCoordinates< PointSetType >( meshSize ); 
      expectedNumPts = meshSize;
      break;
    default:
      std::cerr << "Exit wrong arguments\n";
      std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]\n";
      return EXIT_FAILURE;
      break;
  }
  
  PointSetType::Pointer myPointSet = PointSetType::New();
  for( int i = 0; i < expectedNumPts ; i++ )
    {
    myPointSet->SetPoint( i, pts[i] );
    }
  
  MeshType::Pointer myTriangulatedMesh = MeshType::New();
  try
    {
    myTriangulatedMesh = DelaunayTriangulation< PointSetType, MeshType >( myPointSet );
    }
  catch( int e ) 
    {
    std::cerr << "Exception WiT Caught\n";
    return EXIT_FAILURE;
    }

  CheckONextLink< MeshType >( myTriangulatedMesh );
  typedef itk::QuadEdgeMesh< double, 3 > TempMeshType;
  typedef itk::VTKPolyDataWriter< TempMeshType > MeshWriter;
  MeshWriter::Pointer write = MeshWriter::New();
  write->SetFileName("./tempMesh.vtk");
  write->SetInput( myTriangulatedMesh );
  write->Update();
  

  return EXIT_SUCCESS;
}
