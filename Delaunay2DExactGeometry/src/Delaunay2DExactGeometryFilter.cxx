//--------------
// itk code
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
#include "itkQuadEdgeMeshEulerOperatorFlipEdgeFunction.h"
#include "itkDelaunayConformingQuadEdgeMeshFilter.h"

//---------------
// our code
#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"

//---------------
// std
#include <iostream>
#include <limits>

//--------------
// temporary
#include "itkVTKPolyDataWriter.h"

//--------------------------------------------------------------------------------
// Random coordonates generation function
//
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateRandomPointCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType        PointType;
  typedef typename PointType::CoordRepType CoordRepType;
  std::vector< PointType > oPt( iN * iN );
 
  // NOTE ALEX: is this cross platform? 
  srand(time(NULL));

  for( unsigned int i = 0; i < iN; i++ )
    {
    oPt[ i ][0] = static_cast< CoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][1] = static_cast< CoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][2] = static_cast< CoordRepType >( 0. );
    }
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Dummy point deletion
//
template< class TMeshType >
void
DeleteDummyPoints( TMeshType* mesh )
{
  typedef typename TMeshType::QEType MeshQuadEdgeType;

  for( unsigned int i = 0; i < 4; i++ )
    {
    MeshQuadEdgeType* e;
    while( e =  mesh->GetPoint( i ).GetEdge() )
      {
      mesh->LightWeightDeleteEdge( e );	      
      }
    mesh->DeletePoint( i );
    }
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Dummy Mesh creation function
//
// Create a dummy square mesh for triangulation initialisation
// limit argument is to determine the size of the square
//
template< class TMeshType >
void
CreateDummyMesh( typename TMeshType::Pointer   mesh, 
                 typename TMeshType::PixelType limit )
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
  int i( 0 );  
  std::vector< PointType > pts( 4 );
  
  PixelType min = -limit;  
  PixelType max =  limit;
  
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
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
template<
  class TMeshType
>
typename TMeshType::Pointer
RecursiveFlipEdgeTest( typename TMeshType::Pointer         mesh, 
		       typename TMeshType::PointIdentifier point, 
		       typename TMeshType::CellIdentifier  cell )
{
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
  typedef typename MeshType::QEPrimal                     MeshQuadEdgePrimal;

  typedef typename MeshType::CellsContainerIterator       MeshCellsContainerIteratorType;
  typedef typename MeshType::PointsContainerConstIterator MeshPointsContainerConstIterator;
  
  typedef typename MeshCellType::PointIdConstIterator     MeshCellPointIdConstIterator;
  typedef typename MeshCellType::PointIdIterator          MeshCellPointIdIterator;
  typedef typename MeshCellType::CellAutoPointer          MeshCellCellAutoPointer;

  typedef typename MeshQuadEdgeType::DualOriginRefType    DualOriginRefType;
  
  typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType >                                 QEPolygonCellType;
  typedef typename itk::QuadEdgeMeshEulerOperatorFlipEdgeFunction< MeshType, MeshQuadEdgeType > FlipEdgeFunctionType;
  typedef typename itk::VTKPolyDataWriter< MeshType >                                           MeshWriterType;

  MeshCellCellAutoPointer      cellpointer;
  MeshPointType                pointCoord;
  MeshCellPointIdConstIterator cellPointsIterator;     
  MeshPointIdentifier          q[3];                // points index of T
  MeshPointIdentifier          p[3];                // points index of aT
  int r(0), k(0);                                   // current point and oposite point position  
  
  if( !mesh->GetCell( cell, cellpointer ) )
    {
    std::cerr << "ERROR - Could not find the cell given in parameter." << std::endl;
    throw -1;
    }
  if( !mesh->GetPoint( point, &pointCoord ) )
    {
    std::cerr << "Error - Could not find the point given in parameter." << std::endl;
    throw -1;
    }

  mesh->GetCell(   cell, cellpointer );
  mesh->GetPoint( point, &pointCoord );        
  
  cellPointsIterator = cellpointer->PointIdsBegin();
  int i(0);
  while( cellPointsIterator != cellpointer->PointIdsEnd() )
    {
    p[i] = *cellPointsIterator;
    if( p[i] == point)
      {
      r=i;
      }
    ++cellPointsIterator;
    i++;
    }
      
  MeshQuadEdgeType* e = mesh->FindEdge( p[(r+1)%3], p[(r+2)%3] );
  if( e->IsAtBorder() )
    {
    // the edge is a border edge and can not be flip
    return mesh;
    }

  DualOriginRefType adjCell = e->GetRight();
  if( adjCell == cell )
    {
    adjCell = e->GetLeft();
    }
  mesh->GetCell( adjCell, cellpointer );       
  
  cellPointsIterator = cellpointer->PointIdsBegin();
  int j(0);
  while( cellPointsIterator != cellpointer->PointIdsEnd() )
    {
    q[j] = *cellPointsIterator;
    if( q[j] != p[(r+1)%3] && q[j] != p[(r+2)%3] )
      {
      k = j;
      }
    ++cellPointsIterator;
    j++;
    }
  
   // NOTE STEF: the boolean parameter of the test should be remove in a near future  
   if( TestPointInTriangleInMesh< MeshType >( mesh, adjCell, pointCoord, true ) ) 
     {
     typename FlipEdgeFunctionType::Pointer flipedge = FlipEdgeFunctionType::New();
     flipedge->SetInput( mesh );
     e = flipedge->Evaluate( mesh->FindEdge( p[(r+1)%3], p[(r+2)%3] ) );

     RecursiveFlipEdgeTest< MeshType >( mesh, point, e->GetLeft() );
     RecursiveFlipEdgeTest< MeshType >( mesh, point, e->GetRight() );
     }

  return mesh;
}
//--------------------------------------------------------------------------------

template< class TMeshType >
void
AddPoint( TMeshType* mesh, typename TMeshType::CellIdentifier startingCell, typename TMeshType::PointType point )
{
  typedef          TMeshType                                   MeshType;
  typedef typename MeshType::PointType                         MeshPointType;
  typedef typename MeshType::CellType                          MeshCellType;
  typedef typename MeshType::PointIdentifier                   MeshPointIdentifier;
  typedef typename MeshType::CellIdentifier                    MeshCellIdentifier;
  typedef typename MeshType::PointsContainer                   MeshPointsContainer;
  typedef typename MeshType::CellsContainer                    MeshCellsContainer;
  typedef typename MeshType::PointIdList                       MeshPointIdList;
  typedef typename MeshType::QEType                            MeshQuadEdgeType;
  typedef typename MeshType::CellsContainerIterator            MeshCellsContainerIteratorType;
  typedef typename MeshType::PointsContainerConstIterator      MeshPointsContainerConstIterator;
  typedef typename MeshType::CellAutoPointer                   MeshCellAutoPointer;

  typedef typename MeshCellType::PointIdIterator               MeshCellPointIdIterator;

  typedef typename itk::WalkInTriangulationFunction< MeshType >             WalkInTriangulation;
  typedef typename itk::VectorContainer< unsigned int, MeshCellIdentifier > MeshCellIdVectorContainer;
  typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType >             QEPolygonCellType;


  MeshCellIdentifier  myCellIndex;
  MeshPointIdentifier myPointIndex;
  MeshCellAutoPointer myCellPointer;
  MeshPointType       myPoint;

  WalkInTriangulation* walk = WalkInTriangulation::New();
  MeshCellIdVectorContainer* walkCellIdList = MeshCellIdVectorContainer::New();
  try 
    {
    walkCellIdList = walk->Evaluate( mesh, point, startingCell );
    }
  catch( int e ) 
    {
    std::cerr << "Error - Exception caught in the WalkInTriangulation process" << std::endl;
    throw -1;
    }
  // NOTE STEF: Use a cleaner way to get the last value of a VectorContainer
  myCellIndex = ( --walkCellIdList->End() )->Value();
    
  if( mesh->GetCell( myCellIndex, myCellPointer ) )
    { 
    MeshCellPointIdIterator           pointIdIterator;
    std::vector< MeshCellIdentifier > cellPointsIds( 3 );
    std::vector< MeshCellIdentifier > newCellIds( 3 );
    MeshCellAutoPointer               cellpointer;
    QEPolygonCellType                 *poly;
    
    mesh->GetCell( myCellIndex, myCellPointer );
    
    pointIdIterator = myCellPointer->PointIdsBegin();
    cellPointsIds[0] = *pointIdIterator;
    pointIdIterator++;
    cellPointsIds[1] = *pointIdIterator;
    pointIdIterator++;
    cellPointsIds[2] = *pointIdIterator;
       
    mesh->DeleteFace( myCellIndex );     

    myPointIndex = mesh->FindFirstUnusedPointIndex();
    mesh->SetPoint( myPointIndex, myPoint );

    // NOTE STEF: Use AddTriangleFace method instead
    for( unsigned int i = 0; i < 3; i++ )
      {
      newCellIds[i] = mesh->FindFirstUnusedCellIndex();
      poly = new QEPolygonCellType( 3 );
      cellpointer.TakeOwnership( poly );
      cellpointer->SetPointId( 0, cellPointsIds[ (i)   % 3 ] );
      cellpointer->SetPointId( 1, cellPointsIds[ (i+1) % 3 ] );
      cellpointer->SetPointId( 2, myPointIndex               );
      mesh->SetCell( newCellIds[i], cellpointer );
      } 
        
    RecursiveFlipEdgeTest< MeshType >( mesh, myPointIndex, newCellIds[0] );
    RecursiveFlipEdgeTest< MeshType >( mesh, myPointIndex, newCellIds[1] );
    RecursiveFlipEdgeTest< MeshType >( mesh, myPointIndex, newCellIds[2] );
    }
}	



//--------------------------------------------------------------------------------
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
  
  typedef          TMeshType                                   MeshType;
  typedef typename MeshType::PixelType                         MeshPixelType;
  typedef typename MeshType::Pointer                           MeshTypePointer;
  typedef typename MeshType::PointType                         MeshPointType;
  typedef typename MeshType::CellType                          MeshCellType;
  typedef typename MeshType::PointIdentifier                   MeshPointIdentifier;
  typedef typename MeshType::CellIdentifier                    MeshCellIdentifier;
  typedef typename MeshType::PointsContainer                   MeshPointsContainer;
  typedef typename MeshType::CellsContainer                    MeshCellsContainer;
  typedef typename MeshType::PointIdList                       MeshPointIdList;
  typedef typename MeshType::QEType                            MeshQuadEdgeType;
  typedef typename MeshType::CellsContainerIterator            MeshCellsContainerIteratorType;
  typedef typename MeshType::PointsContainerConstIterator      MeshPointsContainerConstIterator;
  
  typedef itk::VTKPolyDataWriter< MeshType > MeshWriterType;

  typedef typename MeshCellType::PointIdConstIterator          MeshCellPointIdConstIterator;
  typedef typename MeshCellType::PointIdIterator               MeshCellPointIdIterator;
  typedef typename MeshCellType::CellAutoPointer               MeshCellCellAutoPointer;
  
  typedef typename itk::WalkInTriangulationFunction< MeshType > WalkInTriangulation;
  typedef typename WalkInTriangulation::Pointer                 WalkInTriangulationPointer;
  typedef          itk::VectorContainer< unsigned int, int >    MeshCellIdVectorContainer;
  typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType > QEPolygonCellType;

  MeshPixelType infinity = pow( 10, 10 ); // NOTE STEF: This is not infinity
  MeshTypePointer myMesh = MeshType::New();
  CreateDummyMesh< MeshType >( myMesh, infinity );
  
  PointSetPointsContainer              *myPoints           = myPointSet->GetPoints();
  PointSetPointsContainerConstIterator pointIterator       = myPoints->Begin();
  MeshCellIdentifier                   myStartingCellIndex = 0;
  
  int o(1);
  while( pointIterator != myPoints->End() ) 
    { 
      
    PointSetPointType       myTempPoint;
    MeshPointType           myPoint;
    MeshPointIdentifier     myPointIndex;
    MeshCellIdentifier      myCellIndex;
    MeshCellCellAutoPointer myCellPointer;
    
    myTempPoint  = pointIterator.Value();
    myPoint[0]   = myTempPoint[0];
    myPoint[1]   = myTempPoint[1];
    myPoint[2]   = myTempPoint[2];
    
    WalkInTriangulationPointer         letsWalk       = WalkInTriangulation::New();
    MeshCellIdVectorContainer::Pointer walkCellIdList = MeshCellIdVectorContainer::New();
    try 
      {
      walkCellIdList = letsWalk->Evaluate( myMesh, myPoint, myStartingCellIndex );
      }
    catch( int e ) 
      {
      std::cerr << "Error - Exception caught in the WalkInTriangulation process" << std::endl;
      throw -1;
      }
    // NOTE STEF: Use a cleaner way to get the last value of a VectorContainer
    myCellIndex = ( --walkCellIdList->End() )->Value();
    
    if( myMesh->GetCell( myCellIndex, myCellPointer ) )
      { 
      MeshCellPointIdIterator           pointIdIterator;
      std::vector< MeshCellIdentifier > cellPointsIds( 3 );
      std::vector< MeshCellIdentifier > newCellIds( 3 );
      MeshCellCellAutoPointer           cellpointer;
      QEPolygonCellType                 *poly;
      
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

      // NOTE STEF: Use AddTriangleFace method instead
      for( unsigned int i = 0; i < 3; i++ )
        {
        newCellIds[i] = myMesh->FindFirstUnusedCellIndex();
        poly = new QEPolygonCellType( 3 );
        cellpointer.TakeOwnership( poly );
        cellpointer->SetPointId( 0, cellPointsIds[ (i)   % 3 ] );
        cellpointer->SetPointId( 1, cellPointsIds[ (i+1) % 3 ] );
        cellpointer->SetPointId( 2, myPointIndex               );
        myMesh->SetCell( newCellIds[i], cellpointer );
        } 
        
      RecursiveFlipEdgeTest< MeshType >( myMesh, myPointIndex, newCellIds[0] );
      RecursiveFlipEdgeTest< MeshType >( myMesh, myPointIndex, newCellIds[1] );
      RecursiveFlipEdgeTest< MeshType >( myMesh, myPointIndex, newCellIds[2] );
      }
     
    MeshQuadEdgeType* myEdge =  myMesh->GetPoint( myPointIndex ).GetEdge();
    myStartingCellIndex = myEdge->GetLeft(); 
    if( !myMesh->GetCell( myEdge->GetLeft(), myCellPointer ) )
      {
      myStartingCellIndex = myEdge->GetRight(); 
      }

    o++;
    ++pointIterator;
    }    
  
  DeleteDummyPoints< MeshType >( myMesh );

  return myMesh;
}
//--------------------------------------------------------------------------------



//--------------------------------------------------------------------------------
// Main function
//
int
main( int argc, char* argv[] )
{
  if( argc != 3 )
    {
    std::cerr << "Exit wrong arguments" << std::endl;
    std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]";
    std::cerr << std::endl;

    return EXIT_FAILURE;
    }
  
  const unsigned int Dimension = 3;
  typedef double     PixelType;

  typedef itk::QuadEdgeMesh< PixelType, Dimension >                       MeshType; 
  typedef itk::PointSet< PixelType, Dimension >                           PointSetType;
  typedef PointSetType::PointType                                         PointType;
  typedef itk::DelaunayConformingQuadEdgeMeshFilter< MeshType, MeshType > ValidityTestType;
  typedef itk::VTKPolyDataWriter< MeshType >                              MeshWriter;
  
  // -------------------------------------------------- 
  // Initialisation

  int type = atoi( argv[1] );
  int meshSize = atoi( argv[2] );
  int expectedNumPts = 0;
  std::vector< PointType > pts;

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
      std::cerr << "Exit wrong arguments" << std::endl;
      std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]";
      std::cerr << std::endl;
      
      return EXIT_FAILURE;
    }
  
  PointSetType::Pointer myPointSet = PointSetType::New();
  for( int i = 0; i < expectedNumPts ; i++ )
    {
    myPointSet->SetPoint( i, pts[i] );
    }

  // -------------------------------------------------- 
  // Delaunay Test
  
  MeshType::Pointer myTriangulatedMesh = MeshType::New();
  try
    {
    myTriangulatedMesh = DelaunayTriangulation< PointSetType, MeshType >( myPointSet );
    }
  catch( int e ) 
    {
    std::cerr << "Exception WiT Caught" << std::endl;
    return EXIT_FAILURE;
    }

  MeshWriter::Pointer write = MeshWriter::New();
  write->SetFileName("./OutputDelaunayMesh.vtk");
  write->SetInput( myTriangulatedMesh );
  write->Update();

  ValidityTestType::Pointer test = ValidityTestType::New();
  test->SetInput( myTriangulatedMesh );
  test->GraftOutput( myTriangulatedMesh );
  test->Update();

  if( test->GetNumberOfEdgeFlips() > 0 )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
//--------------------------------------------------------------------------------
